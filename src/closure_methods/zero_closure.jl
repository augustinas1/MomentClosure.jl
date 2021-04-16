function zero_closure(sys::MomentEquations, binary_vars::Array{Int,1}=Int[])

    closure = OrderedDict()
    closure_exp = OrderedDict()

    if !isempty(binary_vars)
        sys = bernoulli_moment_eqs(sys, binary_vars)
    end

    if typeof(sys) == CentralMomentEquations
        for i in sys.iter_q
            closure[sys.M[i]] = 0
            closure_exp[sys.M[i]] = 0
        end
    else
        μ = copy(sys.μ)
        μ_symbolic = copy(sys.μ)
        raw_to_central = raw_to_central_moments(sys.N, sys.q_order)
        for i in sys.iter_q
            μ[i] = -(raw_to_central[i]-μ[i])
            closure[sys.μ[i]] = simplify(μ[i], expand=true)
            μ[i] = simplify(substitute(μ[i], closure_exp), expand=true)
            closure_exp[sys.μ[i]] = μ[i]
        end
    end

    close_eqs(sys, closure_exp, closure, true)

end
