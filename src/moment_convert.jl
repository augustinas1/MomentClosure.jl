function cumulants_to_central_moments(N::Int, max_order::Int)

    # obtain cumulants up to (m_order)^th order in terms of
    # central moments using formula from Balakrishan et al. (1998)

    K = Dict()
    M_star = Dict()

    iter_all = construct_iter_all(N, max_order)
    unit_vec = filter(x -> sum(x) == 1, iter_all)
    μ = define_μ(N, 1)
    M = define_M(N, max_order)

    M_star[Tuple(zeros(N))] = 1.0
    for i in 1:N
        eᵢ = unit_vec[i]
        K[eᵢ] = μ[eᵢ]
        M_star[eᵢ] = 0.0
    end

    for order in 2:max_order
        # going up through the orders iteratively to build up the expressions (could do recursively?)

        iter_order = filter(x -> sum(x) == order, iter_all)
        for r in iter_order

            ind = findall(x -> x!= 0, r)[end]
            r_sub = r .- unit_vec[ind]
            iter_i = filter(x -> all(x .<= r_sub), iter_all)
            # find the cumulant \kappa_{\bm{r}}}
            suma = 0.0
            for i in iter_i
                factor = 1.0
                for j in 1:N
                    factor *= binomial(r_sub[j], i[j])
                end
                suma += factor*M[r.-i]*M_star[i]
            end
            K[r] = simplify(suma)

            # Find the central moment \M^*_{\bm{r}}
            suma = 0.0
            for i in iter_i
                factor = 1.0
                for j in 1:N
                    factor *= binomial(r_sub[j], i[j])
                end
                suma += factor*(-K[r.-i])*M_star[i]
            end
            suma -= -μ[unit_vec[ind]]*M_star[r_sub]
            M_star[r] = simplify(suma)

        end

    end

    K

end


function cumulants_to_raw_moments(N::Int, max_order::Int)

    # following Smith (1995)

    #if μ == nothing && typeof(sys) == RawMomentEquations
    #    μ = sys.μ
    #elseif μ == nothing && typeof(sys) == CentralMomentEquations
    #    error("must provide vector μ as it is not contained in CentralMomentEquations")
    #end

    μ = define_μ(N, max_order)
    iter_all = construct_iter_all(N, max_order)
    unit_vec = filter(x -> sum(x) == 1, iter_all)

    K = Dict()
    μ_star = Dict()

    μ_star[Tuple(fill(0, N))] = 1.0
    for i in 1:N
        eᵢ = unit_vec[i]
        K[eᵢ] = μ[eᵢ]
        μ_star[eᵢ] = -μ[eᵢ]
    end

    for order in 2:max_order

        iter_order = filter(x -> sum(x) == order, iter_all)
        for r in iter_order

            ind = findall(x -> x!= 0, r)[end]
            r_sub = r .- unit_vec[ind]
            iter_i = filter(x -> all(x .<= r_sub), iter_all)

            suma = 0.0
            for i in iter_i
                factor = 1.0
                for j in 1:N
                    factor *= binomial(r_sub[j], i[j])
                end
                suma += factor*μ[r.-i]*μ_star[i]
            end
            K[r] = simplify(suma)

            suma = 0.0
            for i in iter_i
                factor = 1.0
                for j in 1:N
                    factor *= binomial(r_sub[j], i[j])
                end
                suma += factor*(-K[r.-i])*μ_star[i]
            end
            μ_star[r] = simplify(suma)
        end

    end

    K

end



function raw_to_central_moments(N::Int, order::Int, μ=nothing; bernoulli=false)

    # Return a dictionary of central moments expressed in terms of raw moments
    # example use:
    # 1 raw_to_central = raw_to_central_moments(2, 3)
    # 2 M₁₂ = raw_to_central[(1,2)] = 2 μ₁₀ μ₀₁² - 2μ₁₁μ₀₁ - μ₀₂μ₁₀ + μ₁₂
    # note that μ is an optional argument which can be used to pass
    # arbitrary values/symbols for each raw moment (different from default μᵢ)

    iter_all = construct_iter_all(N, order)
    iter_μ = filter(x -> sum(x) == 1, iter_all)
    M = define_M(N, order)
    if μ == nothing
        μ = define_μ(N, order)
    elseif typeof(μ) != Dict{Any}{Any} || length(μ) != length(iter_all)
        if bernoulli
            iter_all = keys(μ)
        else
            error("passed arguments are inconsistent (μ vs N & order)")
        end
    else
        iter_all = construct_iter_all(N, order)
    end
    raw_to_central = Dict()

    for i in iter_all

        iter_j = Iterators.filter(x -> all(x .<= i), iter_all)
        suma = 0.0
        for j in iter_j
            term = μ[i.-j]
            for (k, e_k) in zip(1:N, iter_μ)
                term *= (-1)^(j[k])*binomial(i[k], j[k])*μ[e_k]^j[k]
            end
            suma += term
        end
        raw_to_central[i] = simplify(suma)
    end

    raw_to_central

end

function central_to_raw_moments(N::Int, order::Int)

    # Return a dictionary of raw moments expressed in terms of central moments
    # example use:
    # 1 central_to_raw = central_to_raw_moments(2, 3)
    # 2 μ₁₂ = central_to_raw[(1,2)] = 2 M₁₁ μ₀₁ + M₀₂μ₁₀ + M₁₂ + μ₁₀ μ₀₁²

    iter_all = construct_iter_all(N, order)
    iter_μ = filter(x -> sum(x) == 1, iter_all)
    M = define_M(N, order)
    μ = define_μ(N, order)

    central_to_raw = Dict()

    for i in iter_all

        iter_j = Iterators.filter(x -> all(x .<= i), iter_all)
        suma = 0.0
        for j in iter_j
            term = M[i.-j]
            for (k, e_k) in zip(1:N, iter_μ)
                term *= binomial(i[k], j[k])*μ[e_k]^j[k]
            end
            suma += term
        end
        central_to_raw[i] = simplify(suma)
    end

    central_to_raw

end
