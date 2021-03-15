function multi_binomial(a, b)
    # take two vectors a and b as arguments and calculate
    # a product of binomials computed from the corresponding vector components
    # C^a_b = binomial(a₁,b₁) * binomial(a₂,b₂) * ... * binomial(a_n, β_n)
    if length(a) != length(b)
        throw(ErrorException("vector lengths are different for some reason"))
    else
        C = prod([binomial(a[i], b[i]) for i in 1:length(a)])
        return C
    end
end

function derivative_matching(sys::MomentEquations)

    closure_exp = OrderedDict()
    closure = OrderedDict() # additional dict to hold not expanded symbolic expressions

    N = sys.N

    # construct the raw moments up to mth order from the solved-for central moments
    if typeof(sys) == CentralMomentEquations
        closed_μ = central_to_raw_moments(N, sys.m_order)
        μ = central_to_raw_moments(N, sys.q_order)
    else
        μ = sys.μ
        closed_μ = copy(μ)
    end
    μ_symbolic = define_μ(sys.N, sys.q_order)

    # note that derivative matching is originally constructed by truncating at order of m_order+1
    # so if q_order > m_order + 1, we have to consider m_order+1, m_order+2, and so on in sequence
    # to build up the truncated raw moment expressions

    iter_k = vcat(sys.iter_1, sys.iter_m)

    for order in sys.m_order+1:sys.q_order

        # iterator through all moments of lower order
        iter_k = vcat(iter_k, filter(x -> sum(x) == order-1, sys.iter_q))

        length_k = length(iter_k)

        # iterator through all moments of the current truncation order
        iter_order = filter(x -> sum(x) == order, sys.iter_q)

        # initialise matrix needed for linear system defined below
        A = Matrix{Float64}(undef, length_k, length_k)

        for m in iter_order

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
            # While there is no proof that it holds generall, we assume so in order
            # to use integer powers instead of float ones for simpler symbolic manipulations.
            # The code will break otherwise, throwing an InexactError attempting to convert Float to Int.
            # If this ever happens we know what to fix. TODO: add exception handling?
            γ = A\b

            closed_μ[m] = prod([closed_μ[iter_k[i]]^Int(γ[i]) for i in 1:length_k])
            closed_μ[m] = simplify(closed_μ[m])
            closure[μ_symbolic[m]] = prod([μ[iter_k[i]]^Int(γ[i]) for i in 1:length_k])
            closure[μ_symbolic[m]] = simplify(closure[μ_symbolic[m]])
        end

    end

    if typeof(sys) == CentralMomentEquations
        # construct the corresponding truncated expressions of higher order
        # central moments from the obtained raw moment expressions
        raw_to_central = raw_to_central_moments(N, sys.q_order, closed_μ)
        central_to_raw = central_to_raw_moments(N, sys.q_order)
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
