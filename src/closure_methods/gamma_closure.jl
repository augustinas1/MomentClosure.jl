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


function gamma_closure(sys::MomentEquations, binary_vars::Array{Int,1}=Int[])

    closure = OrderedDict() # leaving symbolic higher order terms
    closure_exp = OrderedDict() # expanding out all higher order terms fully

    N = sys.N
    if sys isa CentralMomentEquations
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
            symbolic_sub[αs[i, j]] = M[eᵢ .+ eⱼ] * μ[eᵢ] * μ[eⱼ] / M[2 .* eᵢ] / M[2 .* eⱼ]
        end
        symbolic_sub[x[i]] = μ[eᵢ]^2 / M[2 .* eᵢ]
        β[i] = M[2 .* eᵢ] / μ[eᵢ]
    end

    # construct the raw moments that follow the multivariate gamma distribution
    iters = Vector(undef, N)

    # TupleTools.sort comes in handy here
    unique_iter_q = unique(sort(i) for i in sys.iter_q)
    sub = Dict()

    iter_2nd_order = filter(x-> sum(x) == 2, sys.iter_m)
    # line below reproduces Lakatos et al. (2015) results - stronger assumptions
    #unique_iter_q = unique(sort(i) for i in vcat(sys.iter_m, sys.iter_q))

    for iter in unique_iter_q

        # construct iterator through the product of sums in the Eq. for raw moments
        for j in 1:N
            iter_j = Iterators.filter(x -> sum(x)==iter[j], sys.iter_all)
            iters[j] = iter_j
        end

        iter_k = Iterators.product(iters...)
        suma = 0.0

        multif_iter = multi_factorial(iter)
        for k in iter_k
            factor1 = multif_iter / prod([multi_factorial(k_i) for k_i in k])
            factor2 = 1.0
            for a in 1:N
                factor2 *= gamma_factorial(α[a, a], k[a][a])
                for b in a+1:N
                    factor2 *= gamma_factorial(α[a, b], k[a][b]+k[b][a])
                end
            end
            suma += factor1*factor2
        end

        suma = simplify(suma, expand=true)
        suma = substitute(suma, symbolic_sub)
        suma *= prod([β[j]^iter[j] for j in 1:N])
        μ[iter] = simplify(suma, expand=true)
        closure[μ_symbolic[iter]] = μ[iter]

        # Note that joint moments of multivariate distributions are symmetric in the sense
        # that, for example, M₀₃₀ is equivalent in form to M₃₀₀, only in their respective expressions
        # indices of all lower order moments have to be permuted correspondingly
        # let M₃₀₀ = μ₁₀₀+M₂₁₀, then M₀₃₀ = μ₀₁₀+M₁₂₀ (same form but indices were changed as 1->2 & 2->1)
        # Utilising the symmetry we can greatly reduce the computation time as no redundant
        # expression constructions/simplifications need to be performed anymore
        # TODO: implement for all closures?

        perms = collect(multiset_permutations(iter, length(iter)))[2:end]

        for iter_perm in perms

            iter_perm_ind = sortperm(sortperm(iter_perm))
            for i in sys.iter_1
                sub[μ_symbolic[i]] = μ_symbolic[i[iter_perm_ind]]
            end
            for i in iter_2nd_order
                sub[M[i]] = M[i[iter_perm_ind]]
            end
            #maybe faster to do substitution step by step without dict?
            iter_perm = Tuple(iter_perm)
            μ[iter_perm] = substitute(μ[iter], sub)
            closure[μ_symbolic[iter_perm]] = μ[iter_perm]
        end

    end

    # construct the corresponding truncated expressions of higher order
    # central moments from the obtained gamma raw moment expressions
    if sys isa CentralMomentEquations
        raw_to_central = raw_to_central_moments(N, sys.q_order, μ)
        central_to_raw = central_to_raw_moments(N, sys.q_order)
        closure_M = OrderedDict()
        for i in sys.iter_q
            closure_exp[M[i]] = raw_to_central[i]
            closure_M[M[i]] = simplify(closure[μ_symbolic[i]]-(central_to_raw[i]-M[i]))
        end
        closure = closure_M
    else
        raw_to_central = raw_to_central_moments(N, 2)
        M_to_μ = [M[i] => raw_to_central[i] for i in filter(x -> sum(x)==2, sys.iter_all)]

        # variance of bernoulli variable is expressed in terms of its mean
        #for ind in binary_vars
    #        eᵢ = sys.iter_1[ind]
    #        push!(M_to_μ, M[2 .* eᵢ] => μ[eᵢ] - μ[eᵢ]^2)
    #    end

        for i in sys.iter_q
            closure_exp[sys.μ[i]] = substitute(μ[i], M_to_μ)
            expr = substitute(closure[sys.μ[i]], M_to_μ)
            closure[sys.μ[i]] = simplify(expr)
        end
    end

    if isempty(binary_vars)
        return close_eqs(sys, closure_exp, closure, true)
    else
        # have to perform the bernoulli simplify in the end for gamma closure
        # as otherwise some symbolic terms do not cancel out and break it

        redundant_iter, redundant_eqs, iter_sub = bernoulli_reduce(sys, binary_vars)

        closed_eqs = Equation[]
        for (i, eq) in enumerate(sys.odes.eqs)
            if !(i in redundant_eqs)

                closed_rhs = substitute(eq.rhs, closure_exp)
                closed_rhs = expand(closed_rhs)
                closed_rhs = substitute(closed_rhs, iter_sub)
                closed_rhs = simplify(closed_rhs)

                push!(closed_eqs, Equation(eq.lhs, closed_rhs))

            end
        end

        iter_rm = intersect(sys.iter_q, redundant_iter)

        if sys isa CentralMomentEquations
            moments = M
        else
            moments = μ_symbolic
        end

        for i in iter_rm
            delete!(closure, moments[i])
        end

        vars = extract_variables(closed_eqs, N, sys.q_order)
        odes = ODESystem(closed_eqs, sys.odes.iv, vars, sys.odes.ps)

        return ClosedMomentEquations(odes, closure, sys)

    end

end


#=function gamma_closure(sys::MomentEquations, binary_vars::Array{Int,1}=Int[])

    closure = OrderedDict() # leaving symbolic higher order terms
    closure_exp = OrderedDict() # expanding out all higher order terms fully

    N = sys.N
    if sys isa CentralMomentEquations
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
            symbolic_sub[αs[i, j]] = M[eᵢ .+ eⱼ] * μ[eᵢ] * μ[eⱼ] / M[2 .* eᵢ] / M[2 .* eⱼ]
        end
        symbolic_sub[x[i]] = μ[eᵢ]^2 / M[2 .* eᵢ]
        β[i] = M[2 .* eᵢ] / μ[eᵢ]
    end

    # construct the raw moments that follow the multivariate gamma distribution
    iters = Vector(undef, N)

    # TupleTools come in handy here
    perms_ind = collect(permutations(collect(1:N)))
    unique_iter_q = unique(sort(i) for i in sys.iter_q)
    # line below reproduces Lakatos et al. (2015) results - stronger assumptions
    #unique_iter_q = unique(sort(i) for i in vcat(sys.iter_m, sys.iter_q))

    for iter in unique_iter_q

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

        suma = simplify(suma, expand=true)
        suma = substitute(suma, symbolic_sub)
        suma *= prod([β[j]^iter[j] for j in 1:N])
        μ[iter] = simplify(suma, expand=true)
        closure[μ_symbolic[iter]] = μ[iter]

        # Note that joint moments of multivariate distributions are symmetric in the sense
        # that, for example, M₀₃₀ is equivalent in form to M₃₀₀, only in their respective expressions
        # indices of all lower order moments have to be permuted correspondingly
        # let M₃₀₀ = μ₁₀₀+M₂₁₀, then M₀₃₀ = μ₀₁₀+M₁₂₀ (same form but indices were changed as 1->2 & 2->1)
        # Utilising the symmetry we can greatly reduce the computation time as no redundant
        # expression constructions/simplifications need to be performed anymore
        # TODO: implement for all closures?

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
    if sys isa CentralMomentEquations
        raw_to_central = raw_to_central_moments(N, sys.q_order, μ)
        central_to_raw = central_to_raw_moments(N, sys.q_order)
        closure_M = OrderedDict()
        for i in sys.iter_q
            closure_exp[M[i]] = raw_to_central[i]
            closure_M[M[i]] = simplify(closure[μ_symbolic[i]]-(central_to_raw[i]-M[i]))
        end
        closure = closure_M
    else
        raw_to_central = raw_to_central_moments(N, 2)
        M_to_μ = [M[i] => raw_to_central[i] for i in filter(x -> sum(x)==2, sys.iter_all)]

        # variance of bernoulli variable is expressed in terms of its mean
        #for ind in binary_vars
    #        eᵢ = sys.iter_1[ind]
    #        push!(M_to_μ, M[2 .* eᵢ] => μ[eᵢ] - μ[eᵢ]^2)
    #    end

        for i in sys.iter_q
            closure_exp[sys.μ[i]] = substitute(μ[i], M_to_μ)
            expr = substitute(closure[sys.μ[i]], M_to_μ)
            closure[sys.μ[i]] = simplify(expr)
        end
    end

    if isempty(binary_vars)
        return close_eqs(sys, closure_exp, closure, true)
    else
        # have to perform the bernoulli simplify in the end for gamma closure
        # as otherwise some symbolic terms do not cancel out and break it

        redundant_iter, redundant_eqs, iter_sub = bernoulli_reduce(sys, binary_vars)

        closed_eqs = Equation[]
        for (i, eq) in enumerate(sys.odes.eqs)
            if !(i in redundant_eqs)

                closed_rhs = substitute(eq.rhs, closure_exp)
                closed_rhs = expand(closed_rhs)
                closed_rhs = substitute(closed_rhs, iter_sub)
                closed_rhs = simplify(closed_rhs)

                push!(closed_eqs, Equation(eq.lhs, closed_rhs))

            end
        end

        iter_rm = intersect(sys.iter_q, redundant_iter)

        if sys isa CentralMomentEquations
            moments = M
        else
            moments = μ_symbolic
        end

        for i in iter_rm
            delete!(closure, moments[i])
        end

        vars = extract_variables(closed_eqs, N, sys.q_order)
        odes = ODESystem(closed_eqs, sys.odes.iv, vars, sys.odes.ps)

        return ClosedMomentEquations(odes, closure, sys)

    end

end=#
