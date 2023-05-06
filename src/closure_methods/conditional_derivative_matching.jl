function conditional_derivative_matching(sys::MomentEquations,
                                         binary_vars::Array{Int,1}=Int[])

    if isempty(binary_vars)
        error("conditional closure does not work if there are no binary species")
    end
    N = sys.N
    sys = bernoulli_moment_eqs(sys, binary_vars)

    # define symbolic raw moment expressions
    μ = sys isa CentralMomentEquations ? define_μ(N, sys.q_order) : copy(sys.μ)

    # closure of higher order raw moments without explicit form of truncated moments
    # e.g. μ₁₄ would still be a function of μ₁₃ even though μ₁₃ is also truncated
    closure_μ = OrderedDict()
    # closure of higher order raw moments explicitly expanding closed lower order moments
    # the form we use in the end when solving the ODEs
    closure_μ_exp = OrderedDict()

    ### perform conditional gaussian closure on raw moments μ ###

    # by nonbernoulli we denote moments which cannot be written in the conditional form
    nonbernoulli_iters = filter(x -> sum(x[binary_vars]) == 0, sys.iter_all)

    # iterator through all moments of lower order
    iter_qs = vcat(sys.iter_1, sys.iter_m)
    sub = Dict()

    for order in sys.m_order+1:sys.q_order

        # building the closed moment expressions order by order (due to such hierarchical functional dependency)
        iter_order = filter(x -> sum(x) == order, sys.iter_q)

        for iter in unique(sort(i) for i in iter_order)

            # Two cases: (i) conditional (⟨gp^j⟩=⟨p^j|g=1⟩⟨g⟩) or (ii) basic non bernoulli (⟨p^j⟩)"
            # Here g is a bernoulli variable (e.g gene state), p is any other stochastic variable (e.g no. of proteins)
            if sum(iter[binary_vars]) > 0
                # Case (i): ⟨gp^j⟩=⟨p^j|g=1⟩⟨g⟩
                # boils down to manipulation on indices:
                # 1. ignore the bernoulli variable state and write ⟨gp^j⟩ as μ₀ⱼ
                # 2. perform derivative matching on μ₀ⱼ (assuming cumulant κ₀ⱼ = 0) and multiply by ⟨g⟩ = μ₁₀
                # 3. add back bernoulli variable indices, for example:
                #   moment ⟨gp^3⟩ = ⟨p^3|g=1⟩⟨g⟩
                # i)  ⟨gp^3⟩ = μ₁₃ -> μ₀₃
                # ii) μ₀₃ = 3μ₀₂²*μ₀₁-2μ₀₁^3
                # iii) 3μ₀₂²*μ₀₁-2μ₀₁^3 -> 3μ₁₂²*μ₁₁-2μ₁₁^3

                # step 1
                bernoulli_iter = fill(0, N)
                for ind in binary_vars
                    if iter[ind] > 0
                        bernoulli_iter[ind] = 1
                    end
                end
                bernoulli_iter = Tuple(bernoulli_iter)
                r = iter.-bernoulli_iter

                # step 2

                iter_k = filter(x -> sum(x) < sum(r), iter_qs)
                length_k = length(iter_k)

                A = Matrix{Float64}(undef, length_k, length_k)
                b = [multi_binomial(r, iter_k[i]) for i in 1:length_k]
                for i in 1:length_k
                    for j in 1:length_k
                        A[i, j] = multi_binomial(iter_k[j], iter_k[i])
                    end
                end

                γ = A\b
                conditional_μ = prod([μ[iter_k[i]]^Int(γ[i]) for i in 1:length_k])
                conditional_μ = μ[bernoulli_iter]*conditional_μ
                conditional_μ = simplify(conditional_μ, expand=true)

                # step 3
                iter_conditional = filter(x -> 0 < sum(x) <= sum(r), nonbernoulli_iters)
                #conditional_sub = Dict([Pair(μ[iter], μ[iter.+bernoulli_iter]/μ[bernoulli_iter])
                #                       for iter in iter_conditional])
                conditional_sub = Dict([Pair(μ[iter], μ[iter.+bernoulli_iter]*μ[bernoulli_iter]^-1)
                                       for iter in iter_conditional])
                conditional_μ = substitute(conditional_μ, conditional_sub)
                conditional_μ = simplify(conditional_μ)

                closure_μ[μ[iter]] = conditional_μ
                closure_μ_exp[μ[iter]] = substitute(conditional_μ, closure_μ_exp)
                closure_μ_exp[μ[iter]] = simplify(closure_μ_exp[μ[iter]])
            else
                # Case (ii)

                # Conditional derivative matching has no well-defined approximation for terms
                # such as ⟨p^j⟩ which are independent of the Bernoulli variables
                # We assume that ⟨p^j⟩ are truncated using derivative matching

                iter_k = filter(x -> sum(x) < sum(iter), iter_qs)
                length_k = length(iter_k)

                A = Matrix{Float64}(undef, length_k, length_k)
                b = [multi_binomial(iter, iter_k[i]) for i in 1:length_k]
                for i in 1:length_k
                    for j in 1:length_k
                        A[i, j] = multi_binomial(iter_k[j], iter_k[i])
                    end
                end
                γ = A\b
                moment = prod([μ[iter_k[i]]^Int(γ[i]) for i in 1:length_k])
                moment = simplify(moment)

                closure_μ[μ[iter]] = moment
                closure_μ_exp[μ[iter]] = substitute(moment, closure_μ_exp)
                closure_μ_exp[μ[iter]] = simplify(closure_μ_exp[μ[iter]])

            end

            perms = collect(multiset_permutations(iter, length(iter)))[2:end]

            for iter_perm in perms

                iter_perm_ind = sortperm(sortperm(iter_perm))

                for i in iter_qs
                    sub[μ[i]] = μ[i[iter_perm_ind]]
                end

                iter_perm = Tuple(iter_perm)
                closure_μ[μ[iter_perm]] = substitute(closure_μ[μ[iter]], sub)
                closure_μ_exp[μ[iter_perm]] = substitute(closure_μ_exp[μ[iter]], sub)
            end

        end

        iter_qs = vcat(iter_qs, iter_order)

    end

    if sys isa CentralMomentEquations

        central_to_raw = central_to_raw_moments(N, sys.q_order)
        μ_central = Dict()
        for iter in vcat(sys.iter_m, sys.iter_q)
            μ_central[μ[iter]] = central_to_raw[iter]
        end

        μ_M_exp = define_μ(N, 1)
        for i in sys.iter_m
            μ_M_exp[i] = μ_central[μ[i]]
        end
        μ_M = copy(μ_M_exp)
        for i in sys.iter_q
            μ_M[i] = closure_μ[μ[i]]
            μ_M[i] = substitute(μ_M[i], μ_central)
            μ_M[i] = simplify(μ_M[i])
            μ_M_exp[i] = closure_μ_exp[μ[i]]
            μ_M_exp[i] = substitute(μ_M_exp[i], μ_central)
            μ_M_exp[i] = simplify(μ_M_exp[i])
        end

        closure = OrderedDict()
        closure_exp = OrderedDict()
        # construct the corresponding truncated expressions of higher order
        # central moments from the obtained raw moment expressions
        raw_to_central_exp = raw_to_central_moments(N, sys.q_order, μ_M_exp, bernoulli=true)
        raw_to_central = raw_to_central_moments(N, sys.q_order, μ_M, bernoulli=true)
        for i in sys.iter_q
            closure_exp[sys.M[i]] = expand_mod(raw_to_central_exp[i])
            closure[sys.M[i]] = expand_mod(raw_to_central[i])
        end

    else
        closure_exp = closure_μ_exp
        closure = closure_μ
    end

    close_eqs(sys, closure_exp, closure, false)

end
