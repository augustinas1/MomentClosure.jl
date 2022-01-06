function generate_raw_moment_eqs(drift_eqs::AbstractVector{Equation}, diff::AbstractArray{T}, 
                                 m_order::Int, vars, name, ps, iv) where T <: Union{Sym, Num}
    
    N = length(drift_eqs)
    drift = [e.rhs for e in drift_eqs]
    diff_mat = diff*transpose(diff)
    
    smap = Dict(Pair.(vars, 1:length(vars)))
    drift_factors, drift_powers, drift_order = polynomial_propensities(drift, iv, smap)
    diff_factors, diff_powers, diff_order = polynomial_propensities(value.(diff_mat), iv, smap)
    q_order = m_order + max(drift_order - 1, diff_order - 2)
    
    iter_all = construct_iter_all(N, q_order)
    iter_m = filter(x -> 1 < sum(x) <= m_order, iter_all)
    iter_q = filter(x -> m_order < sum(x) <= q_order, iter_all)
    iter_1 = filter(x -> sum(x) == 1, iter_all)
    
    μ = define_μ(iter_all, iv)
    eqs = Equation[]
    
    for i in 1:N
        poly = sum( drift_factors[i][j] * μ[Tuple(drift_powers[i][j])] for j in 1:length(drift_factors[i]) )
        push!(eqs, Differential(iv)(μ[iter_1[i]]) ~ poly)
    end
    
    for iter in iter_m
        
        term1 = 0
        for i in 1:N
            iszero(iter[i]) ? continue :
            gradf = iter[i]
            gradh = iter .- iter_1[i]
            term1 += gradf*sum( μ[gradh .+ Tuple(drift_powers[i][j])] * drift_factors[i][j] 
                                for j in 1:length(drift_factors[i]) )
        end

        term2 = 0
        for i in 1:N
            for j in 1:N
                ( i == j && isone(iter[i]) ) || (iszero(iter[i]) || iszero(iter[j])) ? continue : 
                hji = iter .- (iter_1[i] .+ iter_1[j])
                hji_factor = i == j ? iter[i]*(iter[i]-1) : iter[i]*iter[j]
                term2 += hji_factor * sum( μ[hji .+ Tuple(diff_powers[i, j][k])] * diff_factors[i, j][k] 
                                                for k in 1:length(diff_factors[i, j]))
            end
        end
        
        push!(eqs, Differential(iv)(μ[iter]) ~ term1 + term2/2) 
    end

    odename = Symbol(name, "_raw_moment_eqs_m", m_order)
    odes = ODESystem(eqs, iv, extract_variables(eqs, μ), ps; name=odename)

    RawMomentEquations(odes, μ, N, m_order, q_order, iter_all, iter_m, iter_q, iter_1)
end

generate_raw_moment_eqs(sys::SDESystem, m_order::Int) = generate_raw_moment_eqs(equations(sys), get_noiseeqs(sys), m_order, 
                                                                                states(sys), nameof(sys), parameters(sys), get_iv(sys))

#=
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

=#