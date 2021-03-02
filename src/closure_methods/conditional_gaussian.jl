function conditional_gaussian_closure(sys::MomentEquations,
                                      binary_vars::Vector=[])

    if isempty(binary_vars)
      error("conditional closure does not work if there are no binary species")
    end
    N = sys.N

    sys = bernoulli_moment_eqs(sys, binary_vars)
    #iter_all = sys.iter_all
    # define symbolic raw moment expressions
    μ = typeof(sys) == CentralMomentEquations ? define_μ(N, sys.q_order) : copy(sys.μ)
    # express cumulants in terms of raw moments
    K = cumulants_to_raw_moments(N, sys.q_order)

    # closure of higher order raw moments without explicit form of truncated moments
    # e.g. μ₁₄ would still be a function of μ₁₃ even though μ₁₃ is also truncated
    closure_μ = Dict()
    # closure of higher order raw moments explicitly expanding closed lower order moments
    # the form we use in the end when solving the ODEs
    closure_μ_exp = Dict()

    # perform conditional gaussian closure on raw moments μ

    # by nonbernoulli we denote moments which cannot be written in the conditional form
    nonbernoulli_iters = filter(x -> sum(x[binary_vars]) == 0, sys.iter_all)
    for order in sys.m_order+1:sys.q_order

        # building the closed moment expressions order by order (due to such hierarchical functional dependency)
        iter_order = filter(x -> sum(x) == order, sys.iter_q)

        for iter in iter_order

            # Two cases: (i) conditional (⟨gp^j⟩=⟨p^j|g=1⟩⟨g⟩) or (ii) basic non bernoulli (⟨p^j⟩)"
            # Here g is a bernoulli variable (e.g gene state), p is any other stochastic variable (e.g no. of proteins)
            if sum(iter[binary_vars]) > 0
                # Case (i): ⟨gp^j⟩=⟨p^j|g=1⟩⟨g⟩
                # boils down to manipulation on indices:
                # 1. ignore the bernoulli variable state and write ⟨gp^j⟩ as μ₀ⱼ
                # 2. perform gaussian closure on μ₀ⱼ (assuming cumulant κ₀ⱼ = 0) and multiply by ⟨g⟩ (μ₁₀)
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
                conditional_μ = simplify(-μ[bernoulli_iter]*expand(K[r] - μ[r]))
                # additional simplification to make expressions lighter

                # step 3
                iter_conditional = filter(x -> 0 < sum(x) <= sum(r), nonbernoulli_iters)
                conditional_sub = Dict([Pair(μ[iter], μ[iter.+bernoulli_iter]/μ[bernoulli_iter])
                                       for iter in iter_conditional])
                conditional_μ = simplify(expand(substitute(conditional_μ, conditional_sub)))

                closure_μ[μ[iter]] = conditional_μ
                closure_μ_exp[μ[iter]] = simplify(substitute(conditional_μ, closure_μ_exp))
                closure_μ_exp[μ[iter]] = simplify(expand(closure_μ_exp[μ[iter]]))

            else
                # Case (ii)

                # Conditional Gaussian closure has no well-defined approximation for terms
                # such as ⟨p^j⟩ which are independent of the Bernoulli variables
                # we assume that the marginal distribution P(p) also follows a Gaussian
                # so that ⟨p^j⟩ is truncated according to normal closure

                moment = simplify(-expand(K[iter]-μ[iter]))

                closure_μ[μ[iter]] = moment
                closure_μ_exp[μ[iter]] = simplify(substitute(moment, closure_μ_exp))
                closure_μ_exp[μ[iter]] = simplify(expand(closure_μ_exp[μ[iter]]))

            end
        end
    end

    if typeof(sys) == CentralMomentEquations

        #central_to_raw = central_to_raw_moments(sys, sys.m_order)
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
        #display(μ_M_exp)
        closure = Dict()
        closure_exp = Dict()
        # construct the corresponding truncated expressions of higher order
        # central moments from the obtained raw moment expressions
        raw_to_central_exp = raw_to_central_moments(N, sys.q_order, μ_M_exp, bernoulli=true)
        raw_to_central = raw_to_central_moments(N, sys.q_order, μ_M, bernoulli=true)
        for i in sys.iter_q
            closure_exp[sys.M[i]] = simplify(raw_to_central_exp[i])
            closure[sys.M[i]] = simplify(raw_to_central[i])
        end

    else
        closure_exp = closure_μ_exp
        closure = closure_μ
    end

    close_eqs(sys, closure_exp, closure)

end
