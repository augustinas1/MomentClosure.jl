function gen_iter(n, d)
    # based on https://twitter.com/evalparse/status/1107964924024635392
    iter = NTuple{n, Int}[]
    for x in partitions(d + n, n)
        x = x .- 1
        if all(x .<= d)
            append!(iter, Tuple.(multiset_permutations(x, n)))
        end
    end
    iter
end

function construct_iter_all(N::Int, order::Int)
    mapreduce(vcat, 0:order) do d
        Base.sort(gen_iter(N, d), rev=true)
    end
end

# Trim a string of form "(a, b, c, d, ...)" to "abcd..."
trim_key(expr) = filter(x -> !(isspace(x) || x == ')' || x== '(' || x==','), string(expr))

# Expand a symbolic expression (no binomial expansion)
expansion_rule_mod = @acrule ~x * +(~~ys) => sum(map(y-> ~x * y, ~~ys)) # apply distribution law
expand_mod = Fixpoint(Prewalk(PassThrough(expansion_rule_mod))) # distributes terms until no longer possible
flatten_rule_mod = @rule(~x::isnotflat(+) => flatten_term(+, ~x)) #
flatten_mod = Fixpoint(PassThrough(flatten_rule_mod)) #
expand_expr = Fixpoint(PassThrough(Chain([expand_mod, flatten_mod]))) # apply flatten and distribution until no longer possible
expand_div = PassThrough(@acrule( +(~~xs) / ~a => sum(map(x -> x / ~a, ~~xs))))

function define_μ(iter::AbstractVector, iv::BasicSymbolic)

    indices = map(trim_key, iter)

    μs = OrderedDict{eltype(iter), Any}()
    for (i, idx) in enumerate(iter)
        if sum(idx) == 0
            μs[idx] = 1
        else
            sym_name = Symbol('μ', join(map_subscripts(indices[i])))
            sym_raw = Sym{FnType{Tuple{Any}, Real}}(sym_name)
            term_raw = Term{Real}(sym_raw, [iv])
            μs[idx] = setmetadata(term_raw, Symbolics.VariableSource,
                                  (:momentclosure, sym_name))
        end
    end

    μs

end

define_μ(N::Int, order::Int, iv::BasicSymbolic) = define_μ(construct_iter_all(N, order), iv)

function define_μ(N::Int, order::Int, iter = construct_iter_all(N, order))
    @parameters t
    return define_μ(iter, value(t))
end


function define_M(iter::AbstractVector, iv::BasicSymbolic)

    indices = map(trim_key, iter)

    Ms = OrderedDict{eltype(iter), Any}()
    for (i, idx) in enumerate(iter)
        if sum(idx) == 0
            Ms[idx] = 1
        elseif sum(idx) == 1
            Ms[idx] = 0
        else
            sym_name = Symbol('M', join(map_subscripts(indices[i])))
            sym_raw = Sym{FnType{Tuple{Any}, Real}}(sym_name)
            term_raw = Term{Real}(sym_raw, [iv])
            Ms[idx] = setmetadata(term_raw, Symbolics.VariableSource,
                                  (:momentclosure, sym_name))
        end
    end

    Ms

end

define_M(N::Int, order::Int, iv::BasicSymbolic) = define_M(construct_iter_all(N, order), iv)

function define_M(N::Int, order::Int, iter = construct_iter_all(N, order))
    @parameters t
    return define_M(iter, value(t))
end

function extract_variables(eqs::Array{Equation, 1}, μ, M=[])

    vars = vcat(values(μ)..., values(M)...)
    # extract variables from rhs of each equation
    eq_vars = unique(vcat(get_variables.(eqs)...))
    # get_variables changes the metadata so have to be careful here...
    intersect!(vars, eq_vars)
    # need this as get_variables does not extract var from `Differential(t)(var(t))`
    diff_vars = [var_from_nested_derivative(eq.lhs)[1] for eq in eqs]
    # filter out the unique ones
    unique(vcat(diff_vars, vars))
    # the correct ordering *should* be preserved

end


## Set of functions to deconstruct polynomial propensities ##

#=
    Consider a polynomial propensity function: a₁ = x^2 + c₁*y²*x,
    where x(t) and y(t) are the molecule numbers variables and c₁ is reaction parameter
    First we need to split a₁ into separate terms, i.e. x^2 and c₁*y²*x. Then we
    determine the independent multiplication factor in each term (1 and c₁). Finally,
    we obtain the power each variable is raised to in each term (x to power 2 in term 1;
    x to power 1 and y to power 2 in term 2). Having this information we can proceed
    in constructing raw moment equations. IF any propensity function is non-polynomial
    then the function `polynomial_propensities` will throw an error.
=#

isvar(x, vars) = any(isequal(x), vars)

"""
    Check that the given expression does NOT depend on the given variables `vars` (expr is constant wrt. vars)
    A variable here is defined as a function of the independent variable `iv`, e.g,  X(t) is variable, where t
    is the independent variable
"""

isconstant(expr, vars, iv) = !istree(expr) || (!isvar(expr, vars) && all(arg -> isconstant(arg, vars, iv), arguments(expr)))

function split_factor_pow(expr, iv, vars)
    base, exp = arguments(expr)
    exp isa Int || error("Unexpected exponent: $expr")

    factor, powers = split_factor(base, iv, vars)
    factor ^ exp, powers .* exp
end

function split_factor_mul(expr, iv, vars)
    powers = zeros(Int, length(vars))
    factor = 1

    for arg in arguments(expr)
        factor_arg, power_arg = split_factor(arg, iv, vars)
        factor *= factor_arg
        powers .+= power_arg
    end

    factor, powers
end

function split_factor_div(expr, iv, vars)
    num, denom = arguments(expr)
    isconstant(denom, vars, iv) || error("The denominator $denom in propensity $expr is not constant.")

    factor, powers = split_factor(num, iv, vars)

    factor / denom, powers
end

function split_factor(expr, iv, vars)
    
    if ispow(expr)
        split_factor_pow(expr, iv, vars)
    elseif ismul(expr)
        split_factor_mul(expr, iv, vars)
    elseif isdiv(expr)
        split_factor_div(expr, iv, vars)
    elseif isconstant(expr, vars, iv)
        expr, zeros(length(vars))
    elseif isvar(expr, vars)
        1, map(isequal(expr), vars)
    else
        error("Expression $expr could not be parsed correctly!")
    end
end


function polynomial_propensity_div(expr, iv, vars)
    num, denom = arguments(expr)
    isconstant(denom, vars, iv) || error("The denominator $denom in propensity $expr is not constant.")
    factors, powers = polynomial_propensity(num, iv, vars)
    factors ./denom, powers
end

function polynomial_propensity_add(expr, iv, vars)
    factors = []
    powers = Vector{Int}[]

    for term in arguments(expr)
        factor_term, power_term = try
            split_factor(term, iv, vars)
        catch e
            error("Propensity function $term is non-polynomial? \n" * string(e))
        end
        push!(factors, factor_term)
        push!(powers, power_term)
    end

    factors, powers
end

function polynomial_propensity(expr, iv, vars)
    if isdiv(expr)
        polynomial_propensity_div(expr, iv, vars)
    elseif isadd(expr)
        polynomial_propensity_add(expr, iv, vars)
    else
        factor, powers = split_factor(expr, iv, vars)
        [factor], [powers]
    end
end

function polynomial_propensities(arr::AbstractArray, iv::BasicSymbolic, smap::AbstractDict)
    vars = [ x for (x,_) in Base.sort(collect(smap), by=x->x[2]) ]

    all_factors = Array{Vector}(undef, size(arr))
    all_powers = Array{Vector{Vector{Int}}}(undef, size(arr))

    for (rind, expr) in enumerate(arr)
        expr = expand(expr)
        all_factors[rind], all_powers[rind] = polynomial_propensity(expr, iv, vars)
    end

    max_power = maximum(sum.(vcat(all_powers...)))
    all_factors, all_powers, max_power
end
