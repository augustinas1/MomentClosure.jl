function zero_closure(sys::MomentEquations, binary_vars::Array{Int,1}=Int[])

    closure = OrderedDict()
    closure_exp = OrderedDict()

    if !isempty(binary_vars)
        sys = bernoulli_moment_eqs(sys, binary_vars)
    end

    if sys isa CentralMomentEquations
        for i in sys.iter_q
            closure[sys.M[i]] = 0
            closure_exp[sys.M[i]] = 0
        end
    else
        μ = copy(sys.μ)
        μ_symbolic = copy(sys.μ)
        raw_to_central = raw_to_central_moments(sys.N, sys.q_order)

        unique_iter_q = unique(sort(i) for i in sys.iter_q)
        sub = Dict()

        for i in unique_iter_q
            μ[i] = -(raw_to_central[i]-μ[i])
            closure[sys.μ[i]] = simplify(μ[i], expand=true)
            μ[i] = simplify(substitute(μ[i], closure_exp), expand=true)
            closure_exp[sys.μ[i]] = μ[i]

            perms = collect(multiset_permutations(i, length(i)))[2:end]

            for iter_perm in perms

                iter_perm_ind = sortperm(sortperm(iter_perm))
                for r in sys.iter_all
                    sub[sys.μ[r]] = sys.μ[r[iter_perm_ind]]
                end

                iter_perm = Tuple(iter_perm)
                closure[sys.μ[iter_perm]] = substitute(closure[sys.μ[i]], sub)
                closure_exp[sys.μ[iter_perm]] = substitute(closure_exp[sys.μ[i]], sub)
            end

        end
    end

    close_eqs(sys, closure_exp, closure, true)

end
