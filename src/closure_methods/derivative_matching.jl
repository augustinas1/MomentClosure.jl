function multi_binomial(a, b)
    # take two vectors a and b as arguments and calculate
    # a product of binomials computed from the corresponding vector components
    # C^a_b = binomial(a₁,b₁) * binomial(a₂,b₂) * ... * binomial(a_n, β_n)
    if length(a) != length(b)
        throw(ErrorException("vector lengths are different for some reason"))
    else
        C = prod([binomial(a[i], b[i]) for i in eachindex(a)])
        return C
    end
end

function derivative_matching(sys::MomentEquations, binary_vars::Array{Int,1}=Int[])

    closure_exp = OrderedDict()
    closure = OrderedDict() # additional dict to hold not expanded symbolic expressions

    N = sys.N
    iv = get_iv(sys)

    if isempty(binary_vars)
        isbernoulli = false
    else
        sys = bernoulli_moment_eqs(sys, binary_vars)
        isbernoulli = true
    end

    # construct the raw moments up to mth order from the solved-for central moments
    if sys isa CentralMomentEquations
        M = sys.M
        closed_μ = central_to_raw_moments(N, sys.m_order; iv)
        μ = central_to_raw_moments(N, sys.q_order; iv)
    else
        μ = sys.μ
        closed_μ = copy(μ)
    end
    μ_symbolic = define_μ(N, sys.q_order, iv)

    # note that derivative matching is originally constructed by truncating at order of m_order+1
    # so if q_order > m_order + 1, we have to consider m_order+1, m_order+2, and so on in sequence
    # to build up the truncated raw moment expressions

    # iterator through all moments of lower order
    iter_k = vcat(sys.iter_1, sys.iter_m)
    sub = Dict()

    for order in sys.m_order+1:sys.q_order
        
        length_k = length(iter_k)

        # iterator through all moments of the current truncation order
        iter_order = filter(x -> sum(x) == order, sys.iter_q)

        # initialise matrix needed for linear system defined below
        A = Matrix{Float64}(undef, length_k, length_k)

        for m in unique(sort(i) for i in iter_order)

            # each higher order moment μ_m is a product of lower order moments ∏ᵣ₌₁ᵏ(μ_mᵣ)^γᵣ
            # the coefficients γᵣ (contained in vector γ) are obtained by solving a system
            # of k equations (k - no. of all lower order moments):
            # C^(m)_(mᵤ) = ∑ᵣ₌₁ᵏ γᵣ * C^(mᵣ)_(mᵤ), s = 1, ..., k
            # where C^(a)_(b) is a product of binomials, so that
            # C^(a)_(b) = binomial(a₁, b₁)*binomial(a₂, b₂)*...*binomial(a_n, b_n)
            # hence for each higher order moment we set up a linear system A*γ=b and solve it

            b = [multi_binomial(m, iter_k[i]) for i in 1:length_k]
            for i in 1:length_k
                for j in 1:length_k
                    A[i, j] = multi_binomial(iter_k[j], iter_k[i])
                end
            end

            # The solution γ appears to always be integer-valued (also mentioned in the paper).
            # While there is no proof that it holds generally, we assume so in order
            # to use integer powers instead of float ones for simpler symbolic manipulations.
            # The code will break otherwise, throwing an InexactError attempting to convert Float to Int.
            # If this ever happens we know what to fix. TODO: add exception handling?
            γ = A\b

            closed_μ[m] = prod([closed_μ[iter_k[i]]^Int(γ[i]) for i in 1:length_k])
            closed_μ[m] = simplify(closed_μ[m])
            closure[μ_symbolic[m]] = prod([μ[iter_k[i]]^Int(γ[i]) for i in 1:length_k])
            closure[μ_symbolic[m]] = simplify(closure[μ_symbolic[m]])

            perms = collect(multiset_permutations(m, length(m)))[2:end]

            for iter_perm in perms

                iter_perm_ind = sortperm(sortperm(iter_perm))
                for i in iter_k
                    sub[μ_symbolic[i]] = μ_symbolic[i[iter_perm_ind]]
                end

                if sys isa CentralMomentEquations
                    for i in sys.iter_m
                        sub[M[i]] = M[i[iter_perm_ind]]
                    end
                end

                iter_perm = Tuple(iter_perm)
                closed_μ[iter_perm] = substitute(closed_μ[m], sub)
                closure[μ_symbolic[iter_perm]] = substitute(closure[μ_symbolic[m]], sub)
            end

        end

        iter_k = vcat(iter_k, iter_order)

    end

    if typeof(sys) == CentralMomentEquations
        # construct the corresponding truncated expressions of higher order
        # central moments from the obtained raw moment expressions
        raw_to_central = raw_to_central_moments(N, sys.q_order, μ=closed_μ, bernoulli=isbernoulli, iv=iv)
        central_to_raw = central_to_raw_moments(N, sys.q_order; iv)
        closure_M = OrderedDict()
        for i in sys.iter_q
            closure_exp[sys.M[i]] = raw_to_central[i]
            expr = simplify(central_to_raw[i]-sys.M[i])
            closure_M[sys.M[i]] = simplify(closure[μ_symbolic[i]]-expr)
        end
        closure = closure_M
    else
        for i in sys.iter_q
            closure_exp[sys.μ[i]] = closed_μ[i]
        end
    end

    close_eqs(sys, closure_exp, closure, false)

end
