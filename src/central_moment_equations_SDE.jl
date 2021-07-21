function generate_centeral_moment_eqs(drift_eqs::AbstractVector{Equation}, diff::AbstractArray{T}, m_order::Int, approx_order::Int = 0, ps = [], iv = nothing) where T <: Union{Sym, Num}
    # if q_order == 0 => assume that all drift_eqs are polynomials
    N = length(drift_eqs)
    ps = Set(ps)
    drift = [e.rhs for e in drift_eqs]
    diffusion_matrix = diff*diff'
    if iv === nothing
        for eq in drift_eqs
            if !(eq.lhs isa Number) # assume eq.lhs is either Differential or Number
                iv = ModelingToolkit.iv_from_nested_derivative(eq.lhs)
                break
            end
        end
    end

    vars = OrderedSet()
    for eq in drift_eqs
        ModelingToolkit.collect_vars!(vars, ps, eq.lhs, iv)
    end
    var_vec = Num[v for v in vars]

    if approx_order == 0
       approx_order = max(maximum(degree(drift, ps)), maximum(degree(diffusion_matrix, ps)))
    end
    q_order = m_order + approx_order - 1

    # iterator over all moments from lowest to highest moment order
    iter_all = construct_iter_all(N, q_order)
    # iterator over central moments up to order m
    iter_m = filter(x -> 1 < sum(x) <= m_order, iter_all)
    # iterator over central moments of order greater than m up to q_order
    iter_q = filter(x -> m_order < sum(x) <= q_order, iter_all)
    # iterator over the first order central moments
    iter_1 = filter(x -> sum(x) == 1, iter_all)
    
    # define moment variables (μ standard, M centered)
    μ = define_μ(iter_1, iv)
    μ_vec = [μ[iter] for iter in iter_1]
    M = define_M(iter_all, iv)

    # determine polynomial approximation
    poly_drift = [taylor_expand(f, var_vec, μ_vec, approx_order, false) for f in drift]
    poly_diffusion_matrix = [taylor_expand(f, var_vec, μ_vec, approx_order, false) for f in diffusion_matrix]

    # determine first order moment equations
    mom_eqs = Equation[]
    mono_to_moment = OrderedDict(prod(vars.^iter) => M[iter] for iter in iter_all)
    moms = Num[]
    for iter in iter_1
        mono = sum(vars .* iter)
        push!(moms, μ[iter])
        push!(mom_eqs, Differential(iv)(μ[iter]) ~ poly_subs(sum(poly_drift .* gradient(mono, vars)), mono_to_moment, ps, true))
    end
    for iter in iter_m
        mono = prod(vars .^ iter)
        push!(moms, M[iter])
        poly = sum(poly_drift .* gradient(mono, vars)) + 1/2* sum( (poly_diffusion_matrix*hessian(mono, vars))[i,i] for i in 1:N ) 
        offset = 0
        for i in 1:N
            if iter[i] != 0 
                iter_incr = Tuple(i != j ? iter[j] : iter[j] - 1 for j in 1:N)
                offset += M[iter_incr]*(iter[i] + 1)*mom_eqs[i].rhs
            end
        end
        push!(mom_eqs, Differential(iv)(M[iter]) ~ - offset + poly_subs(poly, mono_to_moment, ps, true))
    end

    return CentralMomentEquations(ODESystem(mom_eqs, iv, moms, ps), μ, M, N, m_order, q_order, iter_all, iter_m, iter_q, iter_1)
end

generate_central_moment_eqs(sys::SDESystem, m_order::Int, approx_order::Int = 0) = generate_centeral_moment_eqs(sys.eqs, sys.noiseeqs, m_order, approx_order, sys.ps, sys.iv)

