struct RawMomentEquations
    """Moment ODEs describing the evolution of raw moments"""
    odes::ODESystem
    """Symbolic variables defining the raw moments"""
    μ::Dict
    """Number of species within the system"""
    N::Int
    """Order of moment equations"""
    m_order::Int
    """Expansion order"""
    q_order::Int
    """Iterator over all index combinations up to order q_order"""
    iter_all::Vector
    """Iterator over all index combinations up to order m_order"""
    iter_m::Vector
    """Iterator over all index combinations of order greater than m_order up to q_order"""
    iter_q::Vector
    """Iterator over index combinations of order 1"""
    iter_1::Vector
end

"""

"""
function generate_raw_moment_eqs(rn::Union{ReactionSystem,ReactionSystemMod},
                                 m_order::Int; combinatoric_ratelaw = true)

    N = numspecies(rn)
    S = get_S_mat(rn)
    a = propensities(rn, combinatoric_ratelaw = combinatoric_ratelaw)

    term_factors, term_powers, poly_order = polynomial_propensities(a, rn)

    q_order = poly_order + m_order - 1

    # iterator over all moments from lowest to highest moment order
    iter_all = construct_iter_all(N, q_order)
    # iterator over raw moments up to order m
    iter_m = filter(x -> 1 < sum(x) <= m_order, iter_all)
    # iterator over raw moments of order rgrater than m up to q_order
    iter_q = filter(x -> m_order < sum(x) <= q_order, iter_all)
    # iterator over the first order moments
    iter_1 = filter(x -> sum(x) == 1, iter_all)

    μ = define_μ(N, q_order)

    dμ = Dict()
    for i in vcat(iter_1, iter_m)
        dμ[i] = 0
        for r = 1:numreactions(rn)
            iter_j = filter(x -> all(x .<= i) && sum(x) <= sum(i) - 1, iter_all)
            for j in iter_j
                factor_j = 1.0
                for k = 1:N
                    factor_j *= expected_coeff(S[k, r], i[k] - j[k]) * binomial(i[k], j[k])
                end
                suma = 0.0
                for k = 1:length(term_factors[r])
                    suma += term_factors[r][k] * μ[j.+Tuple(term_powers[r][k])]
                end
                dμ[i] += factor_j * suma
                dμ[i] = simplify(dμ[i])
            end
        end
    end

    D = Differential(rn.iv)
    eqs = []
    for i in vcat(iter_1, iter_m)
        push!(eqs, D(μ[i]) ~ dμ[i])
    end

    RawMomentEquations(
        ODESystem(eqs),
        μ,
        N,
        m_order,
        q_order,
        iter_all,
        iter_m,
        iter_q,
        iter_1,
    )

end
