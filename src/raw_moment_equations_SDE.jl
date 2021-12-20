function generate_raw_moment_eqs(drift_eqs::AbstractVector{Equation}, diff::AbstractArray{T}, 
                                 m_order::Int, vars, name, ps, iv) where T <: Union{Sym, Num}
    
    N = length(drift_eqs)
    drift = [e.rhs for e in drift_eqs]
    diff_mat = diff*transpose(diff)
    
    q_order = m_order + max(maximum(degree(eq, vars) for eq in drift) - 1, maximum(degree(eq, vars) for eq in diff_mat) - 2)
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
        H = hessian(mono, vars)
        poly = sum(drift .* gradient(mono, vars)) + 1/2*sum(diff_mat[i,:]'*H[:,i] for i in 1:N )
        push!(mom_eqs, Differential(iv)(mono_to_moment[mono]) ~ poly_subs(poly, mono_to_moment, ps, true))
    end

    odename = Symbol(name, "_raw_moment_eqs_m", m_order)
    odes = ODESystem(mom_eqs, iv, moms, ps; name=odename)

    RawMomentEquations(odes, μ, N, m_order, q_order, iter_all, iter_m, iter_q, iter_1)
end

generate_raw_moment_eqs(sys::SDESystem, m_order::Int) = generate_raw_moment_eqs(equations(sys), get_noiseeqs(sys), m_order, 
                                                                                states(sys), nameof(sys), parameters(sys), get_iv(sys))