function normal_closure(sys::MomentEquations, binary_vars::Array{Int,1}=Int[])

    closure = OrderedDict()
    closure_exp = OrderedDict()
    N = sys.N

    if !isempty(binary_vars)
        sys = bernoulli_moment_eqs(sys, binary_vars)
    end

    # build symbolic expressions of cumulants up to q_order in terms of central/raw moments
    if typeof(sys) == CentralMomentEquations
        moments = sys.M
        K = cumulants_to_central_moments(N, sys.q_order)
    else
        moments = sys.μ
        K = cumulants_to_raw_moments(N, sys.q_order)
    end

    #unique_iter_q = unique(sort(i) for i in sys.iter_q)
    iter_qs = sys.iter_m
    sub = Dict()

    # construct the corresponding truncated expressions of higher order central moments
    for order in sys.m_order+1:sys.q_order

        iter_r = filter(x -> sum(x) == order, sys.iter_q)

        for r in unique(sort(i) for i in iter_r)
            # the last term in the symbolic expression of cumulant κᵣ is Mᵣ (μᵣ)
            # therefore, as we here set κᵣ = 0, only simple manipulation is needed
            closed_moment = -(K[r]-moments[r])
            closed_moment = simplify(closed_moment, expand=true)

            closure[moments[r]] = closed_moment
            closure_exp[moments[r]] = substitute(closed_moment, closure_exp)
            closure_exp[moments[r]] = simplify(closure_exp[moments[r]], expand=true)

            perms = collect(multiset_permutations(r, length(r)))[2:end]

            for iter_perm in perms

                iter_perm_ind = sortperm(sortperm(iter_perm))
                for i in sys.iter_1
                    sub[sys.μ[i]] = sys.μ[i[iter_perm_ind]]
                end
                for i in iter_qs
                    sub[moments[i]] = moments[i[iter_perm_ind]]
                end

                iter_perm = Tuple(iter_perm)
                closure[moments[iter_perm]] = substitute(closure[moments[r]], sub)
                closure_exp[moments[iter_perm]] = substitute(closure_exp[moments[r]], sub)
            end

        end

        iter_qs = vcat(iter_qs, iter_r)

    end

    close_eqs(sys, closure_exp, closure, true)

end
