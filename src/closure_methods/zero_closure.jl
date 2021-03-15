function zero_closure(sys::MomentEquations)

    closure = OrderedDict()
    closure_exp = OrderedDict()

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
            closure[sys.μ[i]] = simplify(μ[i], polynorm=true)
            μ[i] = simplify(substitute(μ[i], closure_exp), polynorm=true)
            closure_exp[sys.μ[i]] = μ[i]
        end
    end

    close_eqs(sys, closure_exp, closure, true)

end
