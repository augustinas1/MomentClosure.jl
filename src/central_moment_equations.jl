#= An assembly of functions enabling the generation of moment
   equations for any chemical reaction network with any type of
   (infinitely differentiable) propensities up to arbitrary order =#

struct CentralMomentEquations
    """Moment ODEs describing the evolution of central moments"""
    odes::ODESystem
    # TODO: merge μ and M to `moments` to make it simpler
    """Symbolic variables defining the means"""
    μ::Dict
    """Symbolic variables defining the central moments"""
    M::Dict
    """Number of Species"""
    N::Int
    """Parameter variables"""
    ps::Vector
    """Order of moment ODEs"""
    m_order::Int
    """Order of moment expansion"""
    exp_order::Int
    # TODO: try to get rid of all the iterators in the struct
    """Iterator over all index combinations up to order exp_order"""
    iter_all::Vector
    """Iterator over all index combinations up to order m_order"""
    iter_m::Vector
    """Iterator over all index combinations of order greater than m_order up to exp_order"""
    iter_exp::Vector
    """Iterator over index combinations of order 1"""
    unit_vec::Vector
end

function fact(i)

    #= Calculate a multi-variate factorial of moment vector i,
       i.e., if i=(a, b, c), then fact(i) = a!b!c! =#

    fact = 1
    for j in i
        fact *= factorial(j)
    end
    return fact

end


function generate_central_moment_eqs(rn::Union{ReactionSystem, ReactionSystemMod},
                                     m_order::Int, exp_order::Int=m_order+1;
                                     combinatoric_ratelaw=true)
    #= Generate the moment equations for the reaction network rn
       up to order of moment expansion m and the order of Taylor
       expansion of propensity functions =#

     # TODO: test this bit is working
    if (m_order >= exp_order)
        error("exp_order must be equal or greater than m_order")
    end

    N = numspecies(rn) # no. of molecular species in the network
    R = numreactions(rn) # no. of reactions in the network

    # propensity functions of all reactions in the network
    a = propensities(rn, combinatoric_ratelaw=combinatoric_ratelaw)

    # net stoichiometric matrix
    S = get_S_mat(rn)

    # iterator over all moments from lowest to highest moment order
    iter_all = construct_iter_all(N, exp_order)
    # iterator over central moments up to order m
    iter_m = filter(x -> 1 < sum(x) <= m_order, iter_all)
    # iterator over central moments of order greater than m up to exp_order
    iter_exp = filter(x -> m_order < sum(x) <= exp_order, iter_all)
    # iterator over the first order central moments
    unit_vec = filter(x -> sum(x) == 1, iter_all)

    #= Define the first order raw moments μ and
       central moments Mᵢ as symbolic variables
       using the functionality of ModelingToolkit.jl =#

    @parameters t # need to define time as a parameter

    μ = define_μ(N, 1)
    μ = delete!(μ, Tuple(zeros(N)))

    M = define_M(N, exp_order)

    #= Obtain all derivatives of the propensity functions with respect
    to all molecular species up to order defined by exp_order.
    For example, Da[2][(1, 2)] is equivalent to (∂³/∂n₁∂n₂²) a₂(μ) =#

    Da = []

    # define the number of molecules of each species as a variable (required for the differentiation)
    n = species(rn)

    # messier code here because no trust in Dict preserving insertion order
    #dict_n_to_μ = #Dict(zip(n, μ)
    dict_n_to_μ = Dict()
    for iter in unit_vec
        ind = findall(x -> x==1, iter)[end]
        dict_n_to_μ[n[ind]] = μ[iter]
    end

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
        #suma = Num(0)
        suma = 0
        for j in iter_all
            suma = suma + Da[r][j]*M[j]*1//fact(j)
            suma = simplify(suma)
            # here // gives fractions (otherwise it's a float number)
            # only issue is that it does not latexify properly
        end

        for i in 1:N
            if r == 1
                push!(du, S[i,r]*suma)
                #du[i] = S[i, r]*suma
            else
                du[i] = S[i, r]*suma + du[i]
            end
            du[i] = simplify(du[i])
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
                iter_k = filter(x -> sum(x) <= exp_order-sum(j), iter_all)
                suma = 0.0
                for k in iter_k
                    suma += Da[r][k]*M[j.+k]*1//fact(k)
                end
                dM[i] += factor_j*suma
                dM[i] = simplify(dM[i])
            end
        end
        for j in 1:N
            if i[j] > 0
                dM[i] -= i[j]*du[j]*M[i.-unit_vec[j]]
            end
        end
        dM[i] = simplify(dM[i])
    end


    #@derivatives D'~t
    D = Differential(t)
    eqs = []
    for i in 1:N
        push!(eqs, D(μ[unit_vec[i]]) ~ du[i])
    end
    for i in iter_m
        push!(eqs, D(M[i]) ~ dM[i])
    end

    return CentralMomentEquations(ODESystem(eqs), μ, M, N, rn.ps, m_order, exp_order,
                                  iter_all, iter_m, iter_exp, unit_vec)
end
