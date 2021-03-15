function log_normal_closure(sys::MomentEquations)

    closure = OrderedDict()
    closure_exp = OrderedDict() # here it does not play a role

    N = sys.N
    if sys isa CentralMomentEquations
        M = copy(sys.M)
        μ = central_to_raw_moments(N, sys.m_order)
        μ_symbolic = define_μ(N, sys.q_order)
    else
        μ = copy(sys.μ)
        μ_symbolic = copy(μ)
    end

    Σ = Dict()
    for j in 1:N
        for k in j:N
            eⱼ = sys.iter_1[j]
            eₖ = sys.iter_1[k]
            if sys isa CentralMomentEquations
                Σ[(j,k)] = 1. + M[eⱼ .+ eₖ] / μ[eⱼ] / μ[eₖ]
            else
                Σ[(j,k)] = μ[eⱼ .+ eₖ] / μ[eⱼ] / μ[eₖ]
            end
        end
    end

    for i in sys.iter_q
        term = prod([μ[sys.iter_1[j]]^i[j] for j in 1:N])
        for j in 1:N
            for k in j+1:N
                term *= Σ[j,k]^(i[j]*i[k])
            end
            term *= Σ[j,j]^(i[j]*(i[j]-1)/2)
        end
        μ[i] = simplify(term)
        closure[μ_symbolic[i]] = μ[i]
    end

    if typeof(sys) == CentralMomentEquations
        # construct the corresponding truncated expressions of higher order
        # central moments from the obtained log-normal raw moment expressions
        raw_to_central = raw_to_central_moments(N, sys.q_order, μ)
        central_to_raw = central_to_raw_moments(N, sys.q_order)
        closure_M = OrderedDict()
        for i in sys.iter_q
            closure_exp[M[i]] = raw_to_central[i]
            closure_M[M[i]] = closure[μ_symbolic[i]]-(central_to_raw[i]-M[i])
            closure_M[M[i]] = simplify(closure_M[M[i]])
        end
        closure = closure_M
    else
        closure_exp = closure
    end

    close_eqs(sys, closure_exp, closure, false)

end
