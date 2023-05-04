"""
    generate_raw_moment_eqs(rn::ReactionSystem, m_order::Int;
                            langevin::Bool=false, combinatoric_ratelaws::Bool=true, smap=speciesmap(rn))

Given a [`ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem)
return the [`RawMomentEquations`](@ref) of the system generated up to `m_order`.

Notes:
- The expansion order ``q``, denoted by `q_order` throughout the docs, is automatically
  determined from the given polynomial form of the propensity functions, see the
  [tutorial](@ref main_tutorial) and the [theory section](@ref raw_moment_eqs) for
  more details on how `q_order` is obtained.
- if `langevin=true`, instead of the Chemical Master Equation the Chemical Langevin
  Equation (diffusion approximation) is considered, and the moment equations are 
  constructed from the corresponding SDE formulation.
- `combinatoric_ratelaws=true` uses binomials in calculating the propensity functions
  of a `ReactionSystem`, see the notes for [`ModelingToolkit.jumpratelaw`]
  (https://mtk.sciml.ai/stable/systems/ReactionSystem/#ModelingToolkit.jumpratelaw).
  *Note* that this field is irrelevant using `ReactionSystemMod` as then the
  propensities are defined directly by the user.
- `smap` sets the variable ordering in the moment equations (which index corresponds to which species
  in the reaction network). By default, this is consistent with the internal system ordering
  accessible with [`Catalyst.speciesmap`](https://catalyst.sciml.ai/stable/api/catalyst_api/#Catalyst.speciesmap).
"""
function generate_raw_moment_eqs(rn::ReactionSystem, m_order::Int;
                                 langevin::Bool=false, combinatoric_ratelaws::Bool=true, smap=speciesmap(rn))
    
    iv = get_iv(rn)
    N = numspecies(rn)
    S = get_stoichiometry(rn, smap)
    a = propensities(rn; combinatoric_ratelaws)

    if langevin
        drift = S*a
        diff = Num[S[i,k] * a[k]^(1//2) for i in 1:N, k in eachindex(a)]
        
        return generate_raw_moment_eqs(Equation[Differential(iv)(s) ~ d for (s, d) in zip(species(rn), drift)], 
                                       diff, m_order, states(rn), nameof(rn), parameters(rn), iv)
    
    end
        
    term_factors, term_powers, poly_order = polynomial_propensities(a, iv, smap)

    q_order = poly_order + m_order - 1

    # iterator over all moments from lowest to highest moment order
    iter_all = construct_iter_all(N, q_order)
    # iterator over raw moments up to order m
    iter_m = filter(x -> 1 < sum(x) <= m_order, iter_all)
    # iterator over raw moments of order greater than m up to q_order
    iter_q = filter(x -> m_order < sum(x) <= q_order, iter_all)
    # iterator over the first order moments
    iter_1 = filter(x -> sum(x) == 1, iter_all)

    μ = define_μ(iter_all, iv)

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
            end
        end
        dμ[i] = simplify(dμ[i])
    end

    D = Differential(iv)
    eqs = Equation[]
    for i in vcat(iter_1, iter_m)
        push!(eqs, D(μ[i]) ~ dμ[i])
    end

    vars = extract_variables(eqs, μ)
    odename = Symbol(nameof(rn), "_raw_moment_eqs_m", m_order)
    odes = ODESystem(eqs, iv, vars, get_ps(rn); name=odename)

    RawMomentEquations(
        odes,
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
