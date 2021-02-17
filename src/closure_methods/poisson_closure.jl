function poisson_closure(sys::Union{RawMomentEquations, CentralMomentEquations})

    closure = Dict()
    closure_exp = Dict()
    N = sys.N

    # build symbolic expressions of cumulants up to exp_order in terms of central/raw moments
    if typeof(sys) == CentralMomentEquations
        moments = copy(sys.M)
        K = cumulants_to_central_moments(N, sys.exp_order)
    else
        moments = copy(sys.μ)
        K = cumulants_to_raw_moments(N, sys.exp_order)
    end

    # construct the corresponding truncated expressions of higher order central moments
    for order in sys.m_order+1:sys.exp_order

        iter_r = filter(x -> sum(x) == order, sys.iter_exp)
        for r in iter_r
            # the last term in the symbolic expression of cumulant κᵣ is Mᵣ (μᵣ)
            # therefore, as we here set κᵣ = 0, only simple manipulation is needed

            # poisson closure - set diagonal higher order cumulants to the corresponding
            # mean value and mixed higher order cumulants to zero
            if sum(r) in r #diagonality condition
                eᵣ = sys.unit_vec[findfirst(!iszero, r)]
                closed_moment = sys.μ[eᵣ]-(K[r]-moments[r])
            else
                closed_moment = -(K[r]-moments[r])
            end
            closed_moment = simplify(expand(closed_moment))

            closure[moments[r]] = closed_moment
            closure_exp[moments[r]] = substitute(closed_moment, closure_exp)
            closure_exp[moments[r]] = simplify(expand(closure_exp[moments[r]]))
        end

    end

    close_eqs(sys, closure_exp, closure)

end
