#= An assembly of functions enabling the generation of moment
   equations for any chemical reaction network with any type of
   (infinitely differentiable) propensities up to arbitrary order =#

function fact(i)

    #= Calculate a multi-variate factorial of moment vector i,
       i.e., if i=(a, b, c), then fact(i) = a!b!c! =#

    fact = 1
    for j in i
        fact *= factorial(j)
    end
    return fact

end

"""
    generate_central_moment_eqs(rn::Union{ReactionSystem, ReactionSystemMod},
                                m_order::Int, q_order::Int=0;
                                langevin::Bool=false, combinatoric_ratelaw::Bool=true, smap=speciesmap(rn))

Given a [`ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem)
or [`ReactionSystemMod`](@ref), return the [`CentralMomentEquations`](@ref) of the system generated up to `m_order`.

Notes:
- if `q_order` is not specified by the user, it is assumed that the reaction network
  contains *only* polynomial propensity functions and hence `q_order` is determined
  automatically as in [`generate_raw_moment_eqs`](@ref). However, `q_order` must
  be specified if non-polynomial propensities are included. Note that the expansion
  order ``q`` denotes the highest order of central moments which will be included
  in the ODEs [(due to the Taylor expansion of propensity functions)](@ref central_moment_eqs).
- if `langevin=true`, instead of the Chemical Master Equation the Chemical Langevin
  Equation (diffusion approximation) is considered, and the moment equations are 
  constructed from the corresponding SDE formulation.
- `combinatoric_ratelaw=true` uses binomials in calculating the propensity functions
  of a `ReactionSystem`, see the notes for [`ModelingToolkit.jumpratelaw`]
  (https://mtk.sciml.ai/stable/systems/ReactionSystem/#ModelingToolkit.jumpratelaw).
  *Note* that this field is irrelevant using `ReactionSystemMod` as then the
  propensities are defined directly by the user.
- `smap` sets the variable ordering in the moment equations (which index corresponds to which species
  in the reaction network). By default, this is consistent with the internal system ordering
  accessible with [`speciesmap`](@ref).
"""
function generate_central_moment_eqs(rn::Union{ReactionSystem, ReactionSystemMod},
                                     m_order::Int, q_order::Int=0;
                                     langevin::Bool=false, combinatoric_ratelaw::Bool=true, smap=speciesmap(rn))

    N = numspecies(rn)   # no. of molecular species in the network
    R = numreactions(rn) # no. of reactions in the network
    iv = get_iv(rn)      # independent variable (usually time)
    a = propensities(rn; combinatoric_ratelaw) # propensity functions of all reactions in the network
    S = netstoichmat(rn) # net stoichiometric matrix
    S = reordered_netstoichmat(rn, S, smap)

    if langevin 
        drift = S*a
        diff = Num[S[i,k] * a[k]^(1//2) for i in 1:N, k in eachindex(a)]
        
        return generate_central_moment_eqs(Equation[Differential(iv)(s) ~ d for (s, d) in zip(species(rn), drift)], 
                                           diff, m_order, q_order, states(rn), nameof(rn), parameters(rn), iv)
    
    end

    # quite messy way to check whether all propensity functions are polynomials
    # and extract the moment expansion order automatically (if not set by the user)
    if iszero(q_order)
        try
            _, _, poly_order = polynomial_propensities(a, iv, smap)
            q_order = poly_order + m_order - 1
        catch e
            error("non-polynomial rates (OR A BUG): please specify q_order.\n" * string(e))
        end
    end

    # iterator over all moments from lowest to highest moment order
    iter_all = construct_iter_all(N, q_order)
    # iterator over central moments up to order m
    iter_m = filter(x -> 1 < sum(x) <= m_order, iter_all)
    # iterator over central moments of order greater than m up to q_order
    iter_q = filter(x -> m_order < sum(x) <= q_order, iter_all)
    # iterator over the first order central moments
    iter_1 = filter(x -> sum(x) == 1, iter_all)

    #= Define the first order raw moments μ and
       central moments Mᵢ as symbolic variables
       using the functionality of ModelingToolkit.jl =#

    μ = define_μ(iter_1, iv)
    M = define_M(iter_all, iv)

    #= Obtain all derivatives of the propensity functions with respect
    to all molecular species up to order defined by q_order.
    For example, Da[2][(1, 2)] is equivalent to (∂³/∂n₁∂n₂²) a₂(μ) =#

    Da = []
    # define the number of molecules of each species as a variable (required for the differentiation)
    n = species(rn)
    dict_n_to_μ = Dict(zip(n, values(μ)))

    # save derivatives in a dictionary
    derivs = Dict()
    for r in 1:R
        for k in iter_all
            derivs[k] = a[r]
            for i in 1:N
                for d in 1:k[i]
                    derivs[k] = expand_derivatives(Differential(n[i])(derivs[k]), true) # simplify = true/false
                end
            end
            derivs[k] = substitute(derivs[k], dict_n_to_μ)
        end
        push!(Da, copy(derivs))
    end

    # generate the equations for the first raw moments (means)
    du = []
    for r in 1:R
        suma = 0
        for j in iter_all
            suma = suma + Da[r][j]*M[j]*1//fact(j)
            #suma = simplify(suma)
            # here // gives fractions (otherwise it's a float number)
            # only issue is that it does not latexify properly
        end

        for i in 1:N
            if r == 1
                push!(du, S[i,r]*suma)
            else
                du[i] = S[i, r]*suma + du[i]
            end
            du[i] = expand(du[i])
        end
    end

    # generate the equations for the central moments up to order m

    dM = Dict()
    for i in iter_m
        dM[i] = 0
        for r in 1:R
            iter_j = filter(x -> all(x .<= i) && sum(x) <= sum(i)-1, iter_all)
            for j in iter_j
                factor_j = 1.0
                for k in 1:N
                    factor_j *= expected_coeff(S[k, r], i[k]-j[k]) * binomial(i[k], j[k])
                end
                # m+1 -> MA_order
                iter_k = filter(x -> sum(x) <= q_order-sum(j), iter_all)
                suma = 0.0
                for k in iter_k
                    suma += Da[r][k]*M[j.+k]*1//fact(k)
                end
                dM[i] += factor_j*suma
            end
        end
        for j in 1:N
            if i[j] > 0
                dM[i] -= i[j]*du[j]*M[i.-iter_1[j]]
            end
        end
        dM[i] = expand(dM[i])
    end

    D = Differential(iv)
    eqs = Equation[]
    for i in 1:N
        push!(eqs, D(μ[iter_1[i]]) ~ du[i])
    end
    for i in iter_m
        push!(eqs, D(M[i]) ~ dM[i])
    end

    vars = extract_variables(eqs, μ, M)
    odename = Symbol(nameof(rn), "_central_moment_eqs_m", m_order, "_q", q_order)
    odes = ODESystem(eqs, iv, vars, get_ps(rn); name=odename)

    CentralMomentEquations(
        odes,
        μ,
        M,
        N,
        m_order,
        q_order,
        iter_all,
        iter_m,
        iter_q,
        iter_1
    )
end
