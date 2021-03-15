function poisson_closure(sys::MomentEquations)

    closure = OrderedDict()
    closure_exp = OrderedDict()
    N = sys.N

    # build symbolic expressions of cumulants up to q_order in terms of central/raw moments
    if typeof(sys) == CentralMomentEquations
        moments = copy(sys.M)
        K = cumulants_to_central_moments(N, sys.q_order)
    else
        moments = copy(sys.μ)
        K = cumulants_to_raw_moments(N, sys.q_order)
    end

    # construct the corresponding truncated expressions of higher order central moments
    for order in sys.m_order+1:sys.q_order

        iter_r = filter(x -> sum(x) == order, sys.iter_q)
        for r in iter_r
            # the last term in the symbolic expression of cumulant κᵣ is Mᵣ (μᵣ)
            # therefore, as we here set κᵣ = 0, only simple manipulation is needed

            # poisson closure - set diagonal higher order cumulants to the corresponding
            # mean value and mixed higher order cumulants to zero
            if sum(r) in r #diagonality condition
                eᵣ = sys.iter_1[findfirst(!iszero, r)]
                closed_moment = sys.μ[eᵣ]-(K[r]-moments[r])
            else
                closed_moment = -(K[r]-moments[r])
            end
            closed_moment = simplify(closed_moment, polynorm=true)

            closure[moments[r]] = closed_moment
            closure_exp[moments[r]] = substitute(closed_moment, closure_exp)
            closure_exp[moments[r]] = simplify(closure_exp[moments[r]], polynorm=true)
        end

    end

    close_eqs(sys, closure_exp, closure, true)

end
