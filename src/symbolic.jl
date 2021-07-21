function gen_iter(n, d)
    # based on https://twitter.com/evalparse/status/1107964924024635392
    iter = []
    for x in partitions(d + n, n)
        x = x .- 1
        if all(x .<= d)
            ys = Set(multiset_permutations(x, n))
            for y in ys
                push!(iter, Tuple(y))
            end
        end
    end
    iter
end

function construct_iter_all(N::Int, order::Int)

    #Construct an ordered iterator going over all moments
    # sequentially in terms of order

    iters = []
    for d in 0:order
        x = Base.sort(gen_iter(N, d), rev=true)
        push!(iters, x...)
    end

    iters

end

# Trim a string of form "(a, b, c, d, ...)" to "abcd..."
trim_key(expr) = filter(x -> !(isspace(x) || x == ')' || x== '(' || x==','), string(expr))

# Expand a symbolic expression (no binomial expansion)
expansion_rule_mod = @acrule ~x * +(~~ys) => sum(map(y-> ~x * y, ~~ys))
expand_mod = Fixpoint(Prewalk(PassThrough(expansion_rule_mod)))
flatten_rule_mod = @rule(~x::isnotflat(+) => flatten_term(+, ~x))
flatten_mod = Fixpoint(PassThrough(flatten_rule_mod))
expand_expr = Fixpoint(PassThrough(Chain([expand_mod, flatten_mod])))

function define_μ(iter::AbstractVector, iv::Union{Sym,Num})

    indices = String[]
    for i in iter
        push!(indices, trim_key(i))
    end

    μs = OrderedDict()
    for (i, idx) in enumerate(iter)
        if sum(idx) == 0
            μs[idx] = 1
        else
            μs[idx] = Term{Real}(Sym{FnType{Tuple{Any}, Real}}(Symbol("μ$(join(map_subscripts(indices[i]), ""))")), [iv])
        end
    end

    μs

end

define_μ(N::Int, order::Int, iv::Union{Sym,Num}) = define_μ(construct_iter_all(N, order), iv)
function define_μ(N::Int, order::Int) 
    @parameters t
    return define_μ(construct_iter_all(N, order), t.val)
end


function define_M(iter::AbstractVector,  iv::Union{Sym,Num})

    indices = String[]
    for i in iter
        push!(indices, trim_key(i))
    end

    Ms = OrderedDict()
    for (i, idx) in enumerate(iter)
        if sum(idx) == 0
            Ms[idx] = 1
        elseif sum(idx) == 1
            Ms[idx] = 0
        else
            Ms[idx] = Term{Real}(Sym{FnType{Tuple{Any}, Real}}(Symbol("M$(join(map_subscripts(indices[i]), ""))")), [iv])
        end
    end

    Ms

end

define_M(N::Int, order::Int, iv::Union{Sym,Num}) = define_M(construct_iter_all(N, order), iv)
function define_M(N::Int, order::Int) 
    @parameters t
    return define_M(construct_iter_all(N, order), t.val)
end

function extract_variables(eqs::Array{Equation, 1}, μ, M)
    vars = vcat(values(μ)..., values(M)...)
    # extract variables from rhs of each equation
    eq_vars = unique(vcat(get_variables.(eqs)...))
    # need this as get_variables does not extract var from `Differential(t)(var(t))`
    diff_vars = [var_from_nested_derivative(eq.lhs)[1] for eq in eqs]
    # filter out the unique ones
    eq_vars = unique(vcat(eq_vars..., diff_vars...))
    # this should preserve the correct ordering

    vars = intersect!(vars, eq_vars)
    vars

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

function isconstant(expr::Union{Symbolic, Number}, vars::Vector, iv::Sym) :: Bool

    # Check that the given expression does NOT depend on the given variables `vars` (expr is constant wrt. vars)
    # A variable here is defined as a function of the independent variable `iv`, e.g,  X(t) is variable, where t
    # is the independent variable

    if istree(expr)
        args = arguments(expr)
        if length(args) == 1
            if isequal(args[1], iv)
                !isvar(expr, vars)
            else
                isconstant(args[1], vars, iv)
            end
        else
            all((isconstant(arg, vars, iv) for arg in args))
        end
    else
        return true
    end

end

isvar(x::Symbolic, vars::Vector) = any((isequal(x, var) for var in vars))

function extract_pwr(expr::Symbolic, smap::Dict, vars::Vector, powers::Vector)
    args = arguments(expr)
    if isvar(args[1], vars) && args[2] isa Int && args[2] >= 0
        idx = smap[args[1]]
        if powers[idx] != 0
            error(args[1], " occurring multiple times in a monomial is unexpected. In ", expr)
        else
            powers[idx] = args[2]
        end
    else
        error("Unexpected term: ", expr)
    end
end

function extract_mul(expr::Symbolic, smap::Dict, vars::Vector, iv::Sym)

    powers = zeros(Int, length(smap))
    factors = []

    for arg in arguments(expr)
        if isconstant(arg, vars, iv)
            push!(factors, arg)
        elseif operation(arg) == ^
            extract_pwr(arg, smap, vars, powers)
        elseif isvar(arg, vars)
            powers[smap[arg]] = 1
        else
            error("Unexpected operation: ", operation(arg), " in ", expr)
        end
    end

    factor = isempty(factors) ? 1 : prod(factors)
    powers, factor

end


function polynomial_propensities(a::Vector, rn::Union{ReactionSystem, ReactionSystemMod})

    R = numreactions(rn)
    N = numspecies(rn)
    vars = species(rn)
    smap = speciesmap(rn)
    iv = rn.iv

    all_factors = [[] for i = 1:R]
    all_powers = [[] for i = 1:R]

    for (rind, expr) in enumerate(a)

        expr = expand(expr)

        if isconstant(expr, vars, iv)

            push!(all_factors[rind], expr)
            push!(all_powers[rind], zeros(Int, N))

        elseif isvar(expr, vars)

            push!(all_factors[rind], 1)
            push!(all_powers[rind], map(v -> isequal(expr, v), vars))

        elseif operation(expr) == ^ #Symbolics.Pow

            try
                powers = zeros(Int, N)
                extract_pwr(expr, smap, vars, powers)
                push!(all_factors[rind], 1)
                push!(all_powers[rind], powers)
            catch e
                error("Propensity function ", expr, " is non-polynomial? \n" * string(e))
            end

        elseif operation(expr) == * #Symbolics.Mul

            try
                powers, factor = extract_mul(expr, smap, vars, iv)
                push!(all_powers[rind], powers)
                push!(all_factors[rind], factor)
            catch e
                error("Propensity function ", expr, " is non-polynomial?\n" * string(e))
            end

        elseif operation(expr) == + #Symbolics.Add

            for term in arguments(expr)

                if isconstant(term, vars, iv)
                    push!(all_factors[rind], term)
                    push!(all_powers[rind], zeros(Int, N))
                else
                    try
                        op = operation(term)
                        if op == ^
                            powers = zeros(Int, N)
                            extract_pwr(term, smap, vars, powers)
                            push!(all_factors[rind], 1)
                            push!(all_powers[rind], powers)
                        elseif op == *
                            powers, factor = extract_mul(term, smap, vars, iv)
                            push!(all_factors[rind], factor)
                            push!(all_powers[rind], powers)
                        elseif isvar(term, vars)
                            powers = zeros(Int, N)
                            powers[smap[term]] = 1
                            push!(all_factors[rind], 1)
                            push!(all_powers[rind], powers)
                        else
                            error("Unexpected operation: ", op, " in ", term)
                        end
                    catch e
                        error("Propensity function ", expr, " is non-polynomial?\n" * string(e))
                    end
                end
            end

        else

            error("Expression $expr could not be parsed correctly!")

        end

    end

    max_power = maximum(sum.(vcat(all_powers...)))
    all_factors, all_powers, max_power

end

gradient(f, vars) = [expand_derivatives(Differential(v)(f)) for v in vars]
hessian(f, vars) = [expand_derivatives(Differential(w)(Differential(v)(f))) for v in vars, w in vars]

degree(ex::AbstractArray, ps::Set = Set()) = [degree(e, ps) for e in ex]
degree(ex::Num, ps::Set = Set()) = degree(expand(ex).val, ps)
degree(ex::Symbolics.Pow, ps::Set = Set()) = ex ∈ ps ? 0 : ex.exp 
degree(ex::Symbolics.Mul, ps::Set = Set()) = sum(ex.dict[key] for key in keys(ex.dict) if key ∉ ps)
degree(ex::Symbolics.Add, ps::Set = Set()) = maximum(degree(sub_ex, ps) for sub_ex in arguments(ex))
degree(ex::Sym, ps::Set = Set()) = ex ∈ ps ? 0 : 1
degree(ex::Term, ps::Set = Set()) = ex ∈ ps ? 0 : 1 
degree(ex::Number, ps::Set = Set()) = 0 

poly_subs(ex::Union{Num, Number, Symbolics.Mul, Symbolics.Add}, subs::AbstractDict, ps::AbstractArray, flag::Bool = false) = poly_subs(ex, subs, Set(ps), flag) 
poly_subs(ex::Num, subs::AbstractDict, ps::Set = Set(), flag::Bool = false) = poly_subs(expand(ex).val, subs, ps, flag) 
poly_subs(ex::Number, ::AbstractDict, ps::Set = Set(), ::Bool = false) = ex 
poly_subs(ex::Sym, subs::AbstractDict, ::Set, ::Bool = false) = haskey(subs, ex) ? subs[ex] : ex
poly_subs(ex::Term, subs::AbstractDict, ::Set, ::Bool = false) = haskey(subs, ex) ? subs[ex] : ex
poly_subs(ex::Symbolics.Pow, subs::AbstractDict, ::Set, flag::Bool = false) = substitute(ex, subs)
poly_subs(ex::Symbolics.Add, subs::AbstractDict, ps::Set = Set(), flag::Bool = false) = sum(poly_subs(term, subs, ps, flag) for term in arguments(expand(ex)))

function poly_subs(ex::Symbolics.Mul, subs::AbstractDict, ps = Set(), flag::Bool = false) 
    mono = 1
    coeff = ex.coeff
    if flag 
        for (base, exponent) in pairs(ex.dict)
            if base ∉ ps && haskey(subs, base) # very hacky but does the job ...
                mono *= base^exponent
            else
                coeff *= base^exponent
            end
        end
    else
        for (base, exponent) in pairs(ex.dict)
            if base ∉ ps 
                mono *= base^exponent
            else
                coeff *= base^exponent
            end
        end
    end
    return haskey(subs, mono) ? coeff*subs[mono] : ex
end


function taylor_expand(f::Union{Num, Term, Symbolics.Mul, Symbolics.Add, Symbolics.Pow}, vars::AbstractVector, ref_point::AbstractVector, order::Int, aux_flag::Bool = true)
    @assert order >= 0 "Order of Taylor expansion must be non-negative"
    @assert length(ref_point) == length(vars) "Variables and reference point must have the same dimension" 
    N = length(vars)
    Δ = vars .- aux_flag * ref_point
    sub_ref_point = Dict(vars[i] => ref_point[i] for i in 1:N)
    f_taylor = substitute(f, sub_ref_point)
    iter_all = filter(x -> sum(x) > 0, MomentClosure.construct_iter_all(N, order))
    for iter in iter_all
        f_taylor += substitute(derivative(f, vars, iter), sub_ref_point)*prod(Δ .^ iter)/MomentClosure.fact(iter)
    end
    return f_taylor
end

function derivative(f::Union{Num, Term, Symbolics.Mul, Symbolics.Add, Symbolics.Pow}, vars::AbstractVector, order::Tuple)
    for i in 1:length(vars)
        if order[i] != 0
            f = nth_differential(f, vars[i], order[i])
        end
    end
    return f
end

function nth_differential(f::Union{Num, Term, Symbolics.Mul, Symbolics.Add, Symbolics.Pow}, var::Union{Num, Sym, Term}, n::Int) 
    for i in 1:n
        f = expand_derivatives(Differential(var)(f), true)
    end
    return f
end