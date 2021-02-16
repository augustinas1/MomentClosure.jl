# input is a symbolic expression of polynomial form:
# e.g. x^2 + x*y^2 + 10 + x + y^3 + p
# ONLY expressions of this kind can be expected to work
# our goal is to express ⟨x^2 + xy^2 + 10 + x + y^3 + p⟩, where each moment is
# described by the typical notation we use for higher order raw moments, e.g. ⟨x^2⟩ = M₂₀
# We start by splitting the whole polynomial into terms and continue explicitly fixing each
# term which is a product of numbers, parameters and variables
function polynomial_terms(expr, terms)
    expr = expand(simplify(expr))
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

function factorise_term(expr, factors, pwrs, rn)
    if istree(expr)
        args = arguments(expr)
        if length(args) == 1
            if args[1] == rn.iv
                idx = speciesmap(rn)[expr]
                if pwrs[idx] != 0
                    println("Same variable occurring multiple times is unexpected in: ", expr)
                else
                    pwrs[idx] = 1
                end
            else
                println("Unexpected symbolic expression: ", expr)
            end
        else
            op = operation(expr)
            if op == ^
                if length(args) == 2
                    # this is fragile (any change in SymbolicUtils/ModelingToolkit can break it)
                    if typeof(args[1]) == Term{Real} && typeof(args[2]) == Int64
                        idx = speciesmap(rn)[args[1]]
                        if pwrs[idx] != 0
                            println("Same variable occurring multiple times is unexpected in: ", expr)
                        else
                            pwrs[idx] = args[2]
                        end
                    elseif typeof(args[1]) != Term{Real} &&
                           typeof(args[2]) != Term{Real}
                        push!(factors, expr)
                    else
                        println("Unexpected expression with ^: ", expr)
                    end
                else
                    println("Only expressions x^n (two symbols) are supported ", expr)
                end
            elseif op == *
                for arg in args
                    factorise_term(arg, factors, pwrs, rn)
                end
            else
                println("Operation in a polynomial term is unexpected: ", op)
            end
        end
    else
        push!(factors, expr)
    end

end


struct RawMomentEquations
    """Moment ODEs describing the evolution of raw moments"""
    odes::ODESystem
    """Symbolic variables defining the raw moments"""
    μ::Dict
    """Number of Species"""
    N::Int
    """Parameter variables"""
    ps::Vector
    """Order of moment ODEs"""
    m_order::Int
    """Order of moment expansion"""
    exp_order::Int
    """Iterator over all index combinations up to order exp_order"""
    iter_all::Vector
    """Iterator over all index combinations up to order m_order"""
    iter_m::Vector
    """Iterator over all index combinations of order greater than m_order up to exp_order"""
    iter_exp::Vector
    """Iterator over index combinations of order 1"""
    unit_vec::Vector
end


function generate_raw_moment_eqs(rn::Union{ReactionSystem,ReactionSystemMod},
                                 m_order::Int; combinatoric_ratelaw = true)

    N = numspecies(rn)
    R = numreactions(rn)
    S = get_S_mat(rn)
    a = propensities(rn, combinatoric_ratelaw = combinatoric_ratelaw)

    term_factors = [[] for i = 1:length(a)]
    term_pwrs = [[] for i = 1:length(a)]
    max_pwr = 0
    for (ind, expr) in enumerate(a)
        terms = []
        polynomial_terms(expr, terms)
        for term in terms
            factors = []
            pwrs = zeros(Int64, numspecies(rn))
            try
                factorise_term(term, factors, pwrs, rn)
            catch e
                error("non-polynomial propensity was included?\n" * string(e))
            end
            factor = isempty(factors) ? 1 : prod(factors)
            max_pwr = sum(pwrs) > max_pwr ? sum(pwrs) : max_pwr
            push!(term_factors[ind], factor)
            push!(term_pwrs[ind], pwrs)
        end
    end

    #println("propensity functions are polynomials up to order ", max_pwr)
    exp_order = max_pwr + m_order - 1
    #println("will encounter moments up to order ", exp_order)

    # iterator over all moments from lowest to highest moment order
    iter_all = construct_iter_all(N, exp_order)
    # iterator over raw moments up to order m
    iter_m = filter(x -> 1 < sum(x) <= m_order, iter_all)
    # iterator over raw moments of order rgrater than m up to exp_order
    iter_exp = filter(x -> m_order < sum(x) <= exp_order, iter_all)
    # iterator over the first order moments
    unit_vec = filter(x -> sum(x) == 1, iter_all)

    μ = define_μ(N, exp_order)

    dμ = Dict()
    for i in vcat(unit_vec, iter_m)
        dμ[i] = 0
        for r = 1:R
            iter_j = filter(x -> all(x .<= i) && sum(x) <= sum(i) - 1, iter_all)
            for j in iter_j
                factor_j = 1.0
                for k = 1:N
                    factor_j *= expected_coeff(S[k, r], i[k] - j[k]) * binomial(i[k], j[k])
                end
                suma = 0.0
                for k = 1:length(term_factors[r])
                    suma += term_factors[r][k] * μ[j.+Tuple(term_pwrs[r][k])]
                end
                dμ[i] += factor_j * suma
                dμ[i] = simplify(dμ[i])
            end
        end
    end

    @parameters t
    D = Differential(t)
    eqs = []
    for i in vcat(unit_vec, iter_m)
        push!(eqs, D(μ[i]) ~ dμ[i])
    end

    return RawMomentEquations(
        ODESystem(eqs),
        μ,
        N,
        rn.ps,
        m_order,
        exp_order,
        iter_all,
        iter_m,
        iter_exp,
        unit_vec,
    )

end
