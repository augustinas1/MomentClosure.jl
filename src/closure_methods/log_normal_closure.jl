function log_normal_closure(sys::Union{RawMomentEquations, CentralMomentEquations})

    closure = Dict()
    closure_exp = Dict() # here it does not play a role

    N = sys.N
    if typeof(sys) == CentralMomentEquations
        M = copy(sys.M)
        μ = central_to_raw_moments(N, sys.m_order)
        μ_symbolic = define_μ(N, sys.q_order)
    else
        M = raw_to_central_moments(N, 2)
        μ = copy(sys.μ)
        μ_symbolic = copy(μ)
    end

    Σ = Matrix{Any}(undef, N, N)
    ν = Array{Any}(undef, N)
    for i in 1:N
        # create a key to access variances, i.e., key = \bm{i} = (0, 0, 2, 0)
        eᵢ = sys.unit_vec[i]
        #  construct diagonal elements of Σ (variances) and means ν
        Σ[i, i] = log( 1 + M[eᵢ .* 2] / μ[eᵢ]^2 )
        ν[i] = log(μ[eᵢ])-Σ[i, i]/2
    end

    for i in 1:N
        for j in i+1:N
            # create key to access covariances, i.e., key = \bm{i} = (1, 0, 0, 1)
            eᵢ = sys.unit_vec[i]
            eⱼ = sys.unit_vec[j]
            #key = eᵢ .+ eⱼ
            # construct mixed elements of Σ (covariances)
            #Σ[i, j] = log(1 + M[eᵢ .+ eⱼ] / exp(ν[i] + ν[j] + (Σ[i, i] + Σ[j, j])/2 ))
            Σ[i, j] = log( 1 + M[eᵢ .+ eⱼ] / (μ[eᵢ]*μ[eⱼ]) )
            Σ[j, i] = Σ[i, j]
        end
    end

    # construct the higher order raw moments that follow the log-normal distribution
    for i_tuple in sys.iter_exp # !!! CHANGING THIS to sys.iter_all recovers Lakatos et al. (2015) implementation
        i_vec = collect(i_tuple) # convert to vector (for the linear algebra below)
        μ[i_tuple] = i_vec'*ν + i_vec'*(Σ*i_vec)/2
        μ[i_tuple] = simplify(expand(μ[i_tuple]))
        μ[i_tuple] = simplify(exp(μ[i_tuple]))
        closure[μ_symbolic[i_tuple]] = μ[i_tuple]
    end

    # NOTE that e.g. exp(log(μ₁)) is not simplified to μ₁ automatically by SymbolicUtils...
    # TODO: use log identities to simplify the expressions?

    if typeof(sys) == CentralMomentEquations
        # construct the corresponding truncated expressions of higher order
        # central moments from the obtained log-normal raw moment expressions
        raw_to_central = raw_to_central_moments(N, sys.q_order, μ)
        central_to_raw = central_to_raw_moments(N, sys.q_order)
        closure_M = Dict()
        for i in sys.iter_exp
            closure_exp[M[i]] = raw_to_central[i]
            closure_M[M[i]] = simplify(closure[μ_symbolic[i]]-(central_to_raw[i]-M[i]))
        end
        closure = closure_M
    else
        closure_exp = closure
    end

    close_eqs(sys, closure_exp, closure)

end
