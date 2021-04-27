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

#=
function construct_iter_all(N::Int, order::Int)

    #Construct an ordered iterator going over all moments
    # sequentially in terms of order

    # helper iterator through all possible combinations of
    # N terms each given in range [0, order]
    iter_full = Iterators.product(fill(0:order, N)...)
    iter_all = Iterators.filter(x -> sum(x) == 0, iter_full)
    for i in 1:order
        iter = Iterators.filter(x -> sum(x) == i, iter_full)
        iter_all = Iterators.flatten((iter_all, iter))
    end

    collect(iter_all)

end
=#

# Trim a string of form "(a, b, c, d, ...)" to "abcd..."
trim_key(expr) = filter(x -> !(isspace(x) || x == ')' || x== '(' || x==','), string(expr))

# Expand a symbolic expression (no binomial expansion)
expansion_rule_mod = @acrule ~x * +(~~ys) => sum(map(y-> ~x * y, ~~ys))
expand_mod = Fixpoint(Prewalk(PassThrough(expansion_rule_mod)))
flatten_rule_mod = @rule(~x::isnotflat(+) => flatten_term(+, ~x))
flatten_mod = Fixpoint(PassThrough(flatten_rule_mod))
expand_expr = Fixpoint(PassThrough(Chain([expand_mod, flatten_mod])))

# Note that throughout we use ModelingToolkit to initialise the variables
# but then change their Num type into SymbolicUtils Term{Real} as it's
# easier to handle in symbolic expressions
function define_μ(N::Int, order::Int, iter=construct_iter_all(N, order))

    indices = []
    for i in iter
        push!(indices, trim_key(i))
    end

    @parameters t
    @variables μ[indices](t)

    μs = OrderedDict()
    for (i, idx) in enumerate(iter)
        μs[idx] = sum(idx) == 0 ? 1 : μ[i].val
    end

    μs

end


function define_M(N::Int, order::Int, iter=construct_iter_all(N, order))

    indices = []
    for i in iter
        push!(indices, trim_key(i))
    end

    @parameters t
    @variables M[indices](t)

    Ms = OrderedDict()
    for (i, idx) in enumerate(iter)
        if sum(idx)==0
            Ms[idx] = 1
        elseif sum(idx)==1
            Ms[idx] = 0
        else
            Ms[idx] = M[i].val
        end
    end

    Ms

end


function extract_variables(eqs::Array{Equation, 1}, N::Int, q_order::Int)

    iters = construct_iter_all(N, q_order)
    iter_μ = filter(x -> sum(x) > 0, iters)
    iter_M = filter(x -> sum(x) > 1, iters)

    μs = values(define_μ(N, q_order, iter_μ))
    Ms = values(define_M(N, q_order, iter_M))
    vars = vcat(μs..., Ms...)
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
    NOTE: this code might easily break with a newer ModelingToolkit/SymbolicUtils version
    TODO: try to make this functionality less cumbersome and test more rigorously
    TODO: use AbstractAlgebra MPoly
=#


function polynomial_terms(expr, terms)

    if istree(expr)
        if operation(expr) == +
            for arg in arguments(expr)
                polynomial_terms(arg, terms)
            end
        else
            push!(terms, expr)
        end
    else
        push!(terms, expr)
    end

end

function factorise_term(expr, factors, powers, rn)
    if istree(expr)
        args = arguments(expr)
        if length(args) == 1
            if isequal(args[1], rn.iv)
                idx = speciesmap(rn)[expr]
                if powers[idx] != 0
                    error("Same variable occurring multiple times is unexpected in: ", expr)
                else
                    powers[idx] = 1
                end
            else
                error("Unexpected symbolic expression: ", expr)
            end
        else
            op = operation(expr)
            if op == ^
                if length(args) == 2
                    # this is fragile (any change in SymbolicUtils/ModelingToolkit might break it)
                    if args[1] isa Term{Real} && args[2] isa Int && args[2] >= 0  # x(t)^i where i is Int>0
                        idx = speciesmap(rn)[args[1]]
                        if powers[idx] != 0
                            error("Same variable occurring multiple times is unexpected in: ", expr)
                        else
                            powers[idx] = args[2]
                        end
                    elseif args[1] isa SymbolicUtils.Add{Real,Int64,Dict{Any,Number}}
                        error("binomial expansion is not supported for: ", expr)
                    elseif !isa(args[1], Term{Real}) && !isa(args[2], Term{Real})
                        push!(factors, expr)
                    else
                        error("Unexpected expression with ^: ", expr)
                    end
                else
                    error("Only expressions x^n (two symbols) are supported ", expr)
                end
            elseif op == *
                for arg in args
                    factorise_term(arg, factors, powers, rn)
                end
            else
                error("Operation in a polynomial term is unexpected: ", op)
            end
        end
    else
        push!(factors, expr)
    end

end

function polynomial_propensities(a::Vector, rn::Union{ReactionSystem, ReactionSystemMod})

    R = numreactions(rn)
    term_factors = [[] for i = 1:R]
    term_powers = [[] for i = 1:R]
    max_power = 0
    for (ind, expr) in enumerate(a)
        terms = []
        polynomial_terms(expand(expr), terms)
        for term in terms
            factors = []
            powers = zeros(Int, numspecies(rn))
            try
                factorise_term(term, factors, powers, rn)
            catch e
                error("Non-polynomial propensity was included?\n" * string(e))
            end
            factor = isempty(factors) ? 1 : prod(factors)
            max_power = sum(powers) > max_power ? sum(powers) : max_power
            push!(term_factors[ind], factor)
            push!(term_powers[ind], powers)
        end
    end

    term_factors, term_powers, max_power

end

#=

function is_polynomial(expr::Symbolic, rn::Union{ReactionSystem, ReactionSystemMod})

    # checking whether a given symbolic expression is a polynomial in molecule numbers

    #terms = []
    #polynomial_terms(expand(expr), terms)
    #ks = sort(collect(keys(poly[2])), lt=SymbolicUtils.:<ₑ)

    sym2term, term2sym = SymbolicUtils._dicts()
    SymbolicUtils.labels!((sym2term, term2sym), expr)

    ks = collect(keys(term2sym))
    println(ks)

    for expr in ks[findall(istree.(ks))]

        vars = get_variables(expr)

        if sum(in.(species(rn), Ref(vars))) > 0
            # if length is one we have time-dependent variable (molecule number)
            if length(vars) > 1
                return false
            end
        end

    end

    return true

end

is_polynomial(expr::Real, rn::Union{ReactionSystem, ReactionSystemMod}) = true
=#
