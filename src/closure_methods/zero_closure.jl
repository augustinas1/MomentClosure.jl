function zero_closure(sys::Union{RawMomentEquations, CentralMomentEquations})

    closure = Dict()
    closure_symbolic = Dict()

    if typeof(sys) == CentralMomentEquations
        for i in sys.iter_exp
            closure[sys.M[i]] = 0
            closure_symbolic[sys.M[i]] = 0
        end
    else
        μ = copy(sys.μ)
        μ_symbolic = copy(sys.μ)
        raw_to_central = raw_to_central_moments(sys.N, sys.exp_order)
        for i in sys.iter_exp
            μ[i] = simplify(-(raw_to_central[i]-μ[i]))
            closure_symbolic[sys.μ[i]] = simplify(expand(μ[i]))

            μ[i] = simplify(expand(substitute(μ[i], closure)))
            closure[sys.μ[i]] = μ[i]
        end
    end

    close_eqs(sys, closure, closure_symbolic)

end
