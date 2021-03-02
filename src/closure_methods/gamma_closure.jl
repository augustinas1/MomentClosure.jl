# Rather complicated mathematical formulation of gamma closure makes a computational implementation painful...
# The main issue is that using SymbolicUtils to simplify the cumbersome symbolic expressions can be very slow
# given a higher order truncation and/or a larger reaction network. Q: would SymPy potentially be faster?
# Unfortunately, heavy simplification is a must if we want tractable equations. Most importantly, it appears
# that fully (or to a satisfying degree) simplified expressions for truncated higher order moments do NOT
# depend on negative powers of central moments and hence deterministic ICs can be used (potentially can be
# shown to hold generally if we come up with a simpler analytical expressions for gamma closure). However,
# if the equations are not fully simplified, we will encounter central/raw moments raised to negative powers
# which will lead to NaN values at t=0 if we use deterministic ICs (setting all central moments to zero initially).
# This can be overcome by setting central moments to near-zero initial values,for RawMomentEquations some form of
# symmetry-breaking(?) might potentially work. Nevertheless, sufficient simplification and numerical stability is desired.

function multi_factorial(a)

    # function to compute multivariate factorials of the form
    # \bm{i!} = i_1!i_2!...i_N!

    return prod([factorial(b) for b in a])

end

function gamma_factorial(a, i)

    # function to compute gamma factors of the form
    # (α)_i = α(α+1)...(α+i-1)

    product = 1.0
    for j in 1:i
        product *= a+i-j
    end

    return product

end


function gamma_closure(sys::MomentEquations, clean=true)

    # optional argument clean
    # - clean = true
    #   sacrificing computation time for simpler and numerically stabler expressions,
    #   the equations SHOULD be simplified sufficiently to avoid central moments raised to negative powers
    #   HOWEVER, this may not hold generally and must be check in case-by-case basis
    # - clean = false
    #   sacrificing numerical stability for much faster results which will produce significantly
    #   more cumbersome expressions

    closure = Dict()
    # closure_symbolic saves raw moment closure expressions, potentially of interest when dealing
    # with CentralMomentEquations (identical to closure if RawMomentEquations were passed)
    closure_exp = Dict()

    # gamma closure leads to very convoluted expressions (especially for RawMomentEquations)
    # which SymbolicUtils.simplify() unfortunately struggles to simplify to a satisfying degree

    N = sys.N
    if typeof(sys) == CentralMomentEquations
        M = copy(sys.M)
        μ = central_to_raw_moments(N, sys.m_order)
        μ_symbolic = define_μ(N, sys.q_order)
    else
        # consider covariances explicitly as M (without converting to raw moments)
        # and converting at the end leads to much more tractable expressions
        M = define_M(N, 2)
        μ = copy(sys.μ)
        μ_symbolic = copy(μ)
    end

    # obtain the shape and scale parameters
    # starting with symbolic α initially and substituting the respective
    # central/raw moments after the expressions are simplified
    iter = Iterators.product(fill(1:N, 2)...)
    iter = collect(iter)
    indices = []
    for i in iter
        push!(indices, trim_key(i))
    end
    @parameters α[indices]
    @parameters x[1:N]
    αs = Matrix{Any}(undef, N, N)
    for (ind, (i,j)) in enumerate(iter)
        αs[i, j] = α[ind].val
    end

    for i in 1:N
        for j in i+1:N
            αs[j, i] = αs[i, j]
        end
    end

    for i in 1:N
        αs[i, i] = x[i].val
        for k in 1:N
            if i != k
                αs[i,i] -= αs[i, k]
            end
        end
    end
    α = αs

    symbolic_sub = Dict()
    β = Array{Any}(undef, N)
    for i in 1:N
        eᵢ = sys.iter_1[i]
        for j in i+1:N
            eⱼ = sys.iter_1[j]
            symbolic_sub[αs[i, j]] = M[eᵢ .+ eⱼ] * μ[eᵢ] * μ[eⱼ] / M[2 .* eᵢ] /M[2 .* eⱼ]
        end
        symbolic_sub[x[i]] = μ[eᵢ]^2 / M[2 .* eᵢ]
        β[i] = M[2 .* eᵢ]/μ[eᵢ]
    end

    # construct the raw moments that follow the multivariate gamma distribution
    iters = Vector(undef, N)
    #all_sums = []

    # TupleTools come in handy here
    perms_ind = collect(permutations(collect(1:N)))
    unique_iter_q = unique(sort(i) for i in sys.iter_q)
    # line below reproduces Lakatos
    #unique_iter_q = unique(sort(i) for i in vcat(sys.iter_m, sys.iter_q))

    for iter in unique_iter_q
        #t = @elapsed begin
        #println(iter)
        # construct iterator through the product of sums in the Eq. for raw moments

        for j in 1:N
            iter_j = Iterators.filter(x -> sum(x)==iter[j], sys.iter_all)
            iters[j] = iter_j
        end

        iter_k = Iterators.product(iters...)

        suma = 0.0

        for k in iter_k
            factor1 = multi_factorial(iter) / prod([multi_factorial(k_i) for k_i in k])
            factor2 = 1.0
            for a in 1:N
                factor2 *= gamma_factorial(α[a, a], k[a][a])
                for b in a+1:N
                    factor2 *= gamma_factorial(α[a, b], k[a][b]+k[b][a])
                end
            end
            suma += factor1*factor2
        end
        if clean
            suma = simplify(suma)
            suma = clean_expr(suma)
            suma = substitute(suma, symbolic_sub)
            suma *= prod([β[j]^iter[j] for j in 1:N])
            μ[iter] = simplify(expand_mod(suma))
        else
            suma = substitute(suma, symbolic_sub)
            suma *= prod([β[j]^iter[j] for j in 1:N])
            μ[iter] = simplify(suma)
        end
        # expressions can be simplified significantly but given a larger system
        # that may take too long...
        closure[μ_symbolic[iter]] = μ[iter]

        # Note that joint moments of multivariate distributions are symmetric in the sense
        # that, for example, M₀₃₀ is equivalent in form to M₃₀₀, only in their respective expressions
        # indices of all lower order moments have to be permuted correspondingly
        # let M₃₀₀ = μ₁₀₀+M₂₁₀, then M₀₃₀ = μ₀₁₀+M₁₂₀ (same form but indices were changed as 1->2 & 2->1)
        # Utilising the symmetry we can greatly reduce the computation time as no redundant
        # expression constructions/simplifications need to be performed anymore
        # TODO: implement for all closures
        perms = collect(permutations(iter))
        unique_inds = findfirst.(isequal.(unique(perms)), [perms])
        unique_perms = perms[unique_inds][2:end]
        unique_perms_ind = perms_ind[unique_inds][2:end]
        iter_2nd_order = filter(x-> sum(x) == 2, sys.iter_m)
        for (iter_perm, iter_perm_ind) in zip(unique_perms, unique_perms_ind)
            iter_perm = Tuple(iter_perm)
            sub = Dict()
            for i in sys.iter_1
                sub[μ_symbolic[i]] = μ_symbolic[i[iter_perm_ind]]
            end
            for i in iter_2nd_order
                sub[M[i]] = M[i[iter_perm_ind]]
            end
            μ[iter_perm] = substitute(μ[iter], sub)
            closure[μ_symbolic[iter_perm]] = μ[iter_perm]
        end

    end

    # construct the corresponding truncated expressions of higher order
    # central moments from the obtained gamma raw moment expressions
    if typeof(sys) == CentralMomentEquations
        raw_to_central = raw_to_central_moments(N, sys.q_order, μ)
        central_to_raw = central_to_raw_moments(N, sys.q_order)
        closure_M = Dict()
        for i in sys.iter_q
            # TODO: check if simplify withing raw_to_central is speed bottleneck if clean=false
            closure_exp[M[i]] = raw_to_central[i]
            expr = simplify(central_to_raw[i]-M[i])
            closure_M[M[i]] = simplify(closure[μ_symbolic[i]]-expr)
        end
        closure = closure_M


    else
        raw_to_central = raw_to_central_moments(N, 2)
        M_to_μ = [M[i] => raw_to_central[i] for i in filter(x -> sum(x)==2, sys.iter_all)]
        for i in sys.iter_q
            closure_exp[sys.μ[i]] = substitute(μ[i], M_to_μ)
            expr = substitute(closure[sys.μ[i]], M_to_μ)
            closure[sys.μ[i]] = simplify(expr)
        end
    end

    close_eqs(sys, closure_exp, closure)

end
