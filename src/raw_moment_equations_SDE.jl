function generate_raw_moment_eqs(drift_eqs::AbstractVector{Equation}, diff::AbstractArray{T}, m_order::Int, ps = [], iv = nothing) where T <: Union{Sym, Num}
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
    
    q_order = m_order + max(maximum(degree(drift, ps)) - 1, maximum(degree(diffusion_matrix, ps)) - 2)
    iter_all = construct_iter_all(N, q_order)
    iter_m = filter(x -> 1 < sum(x) <= m_order, iter_all)
    iter_q = filter(x -> m_order < sum(x) <= q_order, iter_all)
    iter_1 = filter(x -> sum(x) == 1, iter_all)
    
    μ = define_μ(iter_all, iv)
    mono_to_moment = OrderedDict(prod(vars.^iter) => μ[iter] for iter in iter_all)
  
    mom_eqs = Equation[]
    moms = []
    for iter in vcat(iter_1, iter_m)
        mono = prod(vars.^iter)
        poly = sum(drift .* gradient(mono, vars)) + 1/2*sum( (diffusion_matrix*hessian(mono, vars))[i,i] for i in 1:N )
        push!(moms, μ[iter])
        push!(mom_eqs, Differential(iv)(mono_to_moment[mono]) ~ poly_subs(poly, mono_to_moment, ps, true))
    end

    return RawMomentEquations(ODESystem(mom_eqs, iv, moms, ps), μ, N, m_order, q_order, iter_all, iter_m, iter_q, iter_1)
end

generate_raw_moment_eqs(sys::SDESystem, m_order::Int) = generate_raw_moment_eqs(sys.eqs, sys.noiseeqs, m_order, sys.ps, sys.iv)