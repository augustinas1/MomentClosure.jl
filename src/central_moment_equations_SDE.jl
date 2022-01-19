function generate_central_moment_eqs(drift_eqs::AbstractVector{Equation}, diff::AbstractArray, 
                                      m_order::Int, q_order::Int, vars, name, ps, iv)
    
    N = length(drift_eqs)
    drift = [e.rhs for e in drift_eqs]
    diff_mat = diff*transpose(diff)
    
    if iszero(q_order)
        try
            smap = Dict(Pair.(vars, 1:length(vars)))
            _, _, drift_order = polynomial_propensities(drift, iv, smap)
            _, _, diff_order = polynomial_propensities(value.(diff_mat), iv, smap)
            q_order = m_order + max(drift_order - 1, diff_order - 2)
        catch e
            error("non-polynomial rates (OR A BUG): please specify q_order.\n" * string(e))
        end
    end

    iter_all = construct_iter_all(N, q_order)
    iter_m = filter(x -> 1 < sum(x) <= m_order, iter_all)
    iter_q = filter(x -> m_order < sum(x) <= q_order, iter_all)
    iter_1 = filter(x -> sum(x) == 1, iter_all)
    
    μ = define_μ(iter_1, iv)
    M = define_M(iter_all, iv)

    dict_vars_to_μ = Dict(zip(vars, values(μ)))
    
    Df = Vector{Dict}(undef, N)
    derivs = Dict()
    for i in 1:N
        for iter in iter_all
            derivs[iter] = drift[i]
            for j in 1:N
                for d in 1:iter[j]
                    derivs[iter] = expand_derivatives(Differential(vars[j])(derivs[iter]), true)
                end
            end
            derivs[iter] = substitute(derivs[iter], dict_vars_to_μ) 
        end
        Df[i] = copy(derivs)
    end

    Dg = Matrix{Dict}(undef, N, N)
    derivs = Dict() 
    for i in 1:N
        for j in 1:N
            for iter in iter_all
                derivs[iter] = diff_mat[i, j]
                for k in 1:N
                    for d in 1:iter[k]
                        derivs[iter] = expand_derivatives(Differential(vars[k])(derivs[iter]), true)
                    end
                end
                derivs[iter] = substitute(derivs[iter], dict_vars_to_μ)
            end
            Dg[i, j] = copy(derivs)
        end
    end
    
    eqs = Equation[]
    for i in 1:N
        poly = sum( Df[i][j]*M[j]//fact(j) for j in iter_all )
        push!(eqs, Differential(iv)(μ[iter_1[i]]) ~ expand(poly))
    end

    for iter in iter_m
    
        iterq2 = filter(x -> sum(x) <= q_order - sum(iter) + 2, iter_all)
        iterq1 = filter(x -> sum(x) <= q_order - sum(iter) + 1, iterq2)
        
        term1 = 0
        for i in 1:N
            iszero(iter[i]) ? continue :
            gradf = iter[i]
            gradh = iter .- iter_1[i]
            term1 += sum( gradf * Df[i][j] * M[gradh .+ j] // fact(j) for j in iterq1 )
            term1 -= sum( gradf * Df[i][j] * M[gradh] * M[j] // fact(j) for j in iter_all ) 
        end

        term2 = 0
        for i in 1:N
            for j in 1:N
                ( i == j && isone(iter[i]) ) || (iszero(iter[i]) || iszero(iter[j])) ? continue : 
                hji = iter .- (iter_1[i] .+ iter_1[j])
                hji_factor = i == j ? iter[i]*(iter[i]-1) : iter[i]*iter[j]

                term2 += hji_factor * sum( Dg[i, j][k] * M[hji .+ k] // fact(k) for k in iterq2 )
            end
        end

        push!(eqs, Differential(iv)(M[iter]) ~ expand(term1 + term2/2))

    end

    odename = Symbol(name, "_central_moment_eqs_m", m_order, "_q", q_order)
    odes = ODESystem(eqs, iv, extract_variables(eqs, μ, M), ps; name=odename)

    CentralMomentEquations(odes, μ, M, N, m_order, q_order, iter_all, iter_m, iter_q, iter_1)
end

"""
    generate_central_moment_eqs(sys::SDESystem, m_order::Int, q_order::Int=0)

Given an [`SDESystem`](https://mtk.sciml.ai/stable/systems/SDESystem/#ModelingToolkit.SDESystem), 
return the [`CentralMomentEquations`](@ref) of the system generated up to `m_order`.
"""
generate_central_moment_eqs(sys::SDESystem, m_order::Int, q_order::Int=0) = generate_central_moment_eqs(equations(sys), get_noiseeqs(sys), m_order, q_order, 
                                                                                                        states(sys), nameof(sys), parameters(sys), get_iv(sys))

#=
function generate_central_moment_eqs(drift_eqs::AbstractVector{Equation}, diff::AbstractArray{T}, 
                                     m_order, approx_order, vars, name, ps, iv) where T <: Union{Sym, Num}
    
    N = length(drift_eqs)
    drift = [e.rhs for e in drift_eqs]
    diffusion_matrix = diff*transpose(diff)

    # if approx_order == 0 => assume that all drift_eqs are polynomials
    if iszero(approx_order)
        approx_order = max(maximum(degree(eq, vars) for eq in drift), maximum(degree(eq, vars) for eq in diffusion_matrix))
        q_order = m_order + max(maximum(degree(eq, vars) for eq in drift) - 1, maximum(degree(eq, vars) for eq in diffusion_matrix) - 2)
    else
        q_order = m_order + approx_order - 1
    end
    
    iter_all = construct_iter_all(N, q_order)
    iter_m = filter(x -> 1 < sum(x) <= m_order, iter_all)
    iter_q = filter(x -> m_order < sum(x) <= q_order, iter_all)
    iter_1 = filter(x -> sum(x) == 1, iter_all)
    
    μ = define_μ(iter_1, iv)
    μ_vec = [value(μ[iter]) for iter in iter_1]
    M = define_M(iter_all, iv)

    # determine polynomial approximation
    # TODO: compare the speed of this + poly_subs to what I have in central_moment_eqs.jl
    poly_drift = [taylor_expand(f, vars, μ_vec, approx_order, false) for f in drift]
    poly_diffusion_matrix = [taylor_expand(f, vars, μ_vec, approx_order, false) for f in diffusion_matrix]

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
        poly = sum(poly_drift .* gradient(mono, vars)) + 1/2 * sum( (poly_diffusion_matrix*hessian(mono, vars))[i,i] for i in 1:N ) 
        offset = 0
        for i in 1:N
            if !iszero(iter[i]) 
                iter_incr = Tuple(i != j ? iter[j] : iter[j] - 1 for j in 1:N)
                offset += M[iter_incr]*(iter[i] + 1)*mom_eqs[i].rhs
            end
        end
 
        push!(mom_eqs, Differential(iv)(M[iter]) ~ - offset + poly_subs(poly, mono_to_moment, Set(ps), true))
    end

    return CentralMomentEquations(ODESystem(mom_eqs, iv, moms, ps, name = Symbol(name, :_RAW_CENTRAL)), μ, M, N, m_order, q_order, iter_all, iter_m, iter_q, iter_1)
end

generate_central_moment_eqs(sys::SDESystem, m_order::Int, q_order::Int=0) = generate_central_moment_eqs(equations(sys), get_noiseeqs(sys), m_order, q_order, 
                                                                                                        states(sys), nameof(sys), parameters(sys), get_iv(sys))
=#