
function sample_raw_moments(sol::EnsembleSolution, order::Int; naive::Bool=true)
    # need to build a t x n matrix for each timepoint, where t is the number of realisations (trajectories)
    # and n is the number of  variables

    # naive algorithm (naive=true), naivemoment(), performs much better in case of small N
    # when N is large, big performance improvement is expected(?) using moment() (naive=false)
    # this implementatation is experimental, need more testing to properly assert any
    # statement about performance + understand why the naive algorithm is still very fast
    # compared to the most trivial algorithm

    # TODO: rewrite this as sample_central_moments right now appearst to be faster
    # TODO: use internal SciMLBase functionality to simplify the code here further?

    alg = naive ? naivemoment : moment

    no_t_pts = length(sol.u[1].t)
    N = length(sol.u[1].u[1])

    iter_all = MomentClosure.construct_iter_all(N, order)
    iter_zero = Tuple(fill(0, N))

    μ = Dict()

    for m in 1:order

        data = Array{Float64}(undef, fill(N, m)..., no_t_pts)

        for t_pt in 1:no_t_pts
            tslice = Array(hcat([s[t_pt] for s in sol]...)')
            data[fill(:, m)..., t_pt] = Array(alg(tslice, m))
        end

        for iter in filter(x-> sum(x) == m, iter_all)
            key = vcat([fill(i, m) for (i,m) in enumerate(iter)]...)
            μ[iter] = data[key..., :]
        end

    end

    μ[iter_zero] = fill(1.0, no_t_pts)

    #typeof(sol) == EnsembleSolution does not work...

    return μ
end


function sample_central_moments(sol::EnsembleSolution, order::Int; naive::Bool=true)

    alg = naive ? naivemoment : moment

    no_t_pts = length(sol.u[1].t)
    N = length(sol.u[1].u[1])

    iter_all = MomentClosure.construct_iter_all(N, order)
    iter_zero = Tuple(fill(0, N))
    iter_order = [filter(x-> sum(x) == m, iter_all) for m in 1:order]

    keys = Dict([iter => vcat([fill(i, n) for (i,n) in enumerate(iter)]...) for iter in iter_all])
    M = Dict([iter => Array{Float64}(undef, no_t_pts) for iter in iter_all])

    for t_pt in 1:no_t_pts
        tslice = Array(hcat([s[t_pt] for s in sol]...)')
        μ = Dict()
        μ[iter_zero] = 1
        for m in 1:order
            data = Array(alg(tslice, m))
            for iter in iter_order[m]
                μ[iter] = data[keys[iter]...]
            end
        end
        M_temp = raw_to_central_moments(N, order, μ)
        for (key, val) in M_temp
            M[key][t_pt] = val
        end
    end

    return M

end


function sample_cumulants(sol::EnsembleSolution, order::Int; naive::Bool=true)

    no_t_pts = length(sol.u[1].t)
    N = length(sol.u[1].u[1])

    iter_all = MomentClosure.construct_iter_all(N, order)
    iter_order = [filter(x-> sum(x) == m, iter_all) for m in 1:order]

    if naive

        κ = Dict()

        for m in 1:order

            data = Array{Float64}(undef, fill(N, m)..., no_t_pts)

            for t_pt in 1:no_t_pts
                tslice = Array(hcat([s[t_pt] for s in sol]...)')
                data[fill(:, m)..., t_pt] = Array(naivecumulant(tslice, m))
            end

            for iter in iter_order[m]
                key = vcat([fill(i, m) for (i,m) in enumerate(iter)]...)
                κ[iter] = data[key..., :]
            end

        end

    else

        iter_all = filter(x -> sum(x) > 0, iter_all)
        keys = Dict([iter => vcat([fill(i, n) for (i,n) in enumerate(iter)]...) for iter in iter_all])
        κ = Dict([iter => Array{Float64}(undef, no_t_pts) for iter in iter_all])

        for t_pt in 1:no_t_pts

            tslice = Array(hcat([s[t_pt] for s in sol]...)')
            κ_tensor = cumulants(tslice, order)

            for m in 1:order

                data = Array(κ_tensor[m])

                for iter in iter_order[m]
                    κ[iter][t_pt] = data[keys[iter]...]
                end

            end

        end

    end

    return κ

end

"""
    deterministic_IC(u₀::Array{T, 1}, eqs::MomentEquations) where T<:Real

Given an array of initial molecule numbers and the corresponding moment equations,
return a mapping of each moment to its initial value under deterministic initial conditions.

# Notes
- The means are set to initial molecule numbers (as they take the values specified in
  `u₀` with probability one). The higher order raw moments are products of the corresponding
  powers of the means whereas the higher order central moments are simply zero.
- The ordering of `u₀` elements must be consistent with the ordering of species
  in the corresponding reaction system (can be checked with the `speciesmap` function).
"""
function deterministic_IC(u₀::Array{T, 1}, eqs::MomentEquations) where T<:Real

    if eqs isa ClosedMomentEquations
        sys = eqs.open_eqs
    else
        sys = eqs
    end

    odes = eqs.odes
    N = sys.N
    if N != length(u₀)
        error("length of the passed IC vector and the number of species in the system are inconsistent")
    end

    μ_map = [sys.μ[iter] => u₀[i] for (i, iter) in enumerate(sys.iter_1)]

    no_states = length(odes.states)
    if typeof(sys) == CentralMomentEquations
        moment_map = [odes.states[i] => 0.0 for i in N+1:no_states]
    else
        reverse_μ = Dict(μ => iter for (iter, μ) in sys.μ)
        moment_map = [odes.states[i] => prod(u₀ .^ reverse_μ[odes.states[i]]) for i in N+1:no_states]
    end

    vcat(μ_map, moment_map)

end

"""
    format_moment_eqs(eqs::MomentEquations)

Given `MomentEquations`, return an array of formatted strings representing
each ODE that can be visualised using [Latexify](https://github.com/korsbo/Latexify.jl).
Although Latexify can be applied directly on [`ModelingToolkit.ODESystem`]
(https://mtk.sciml.ai/stable/systems/ODESystem/) (accessed by `eqs.odes`),
`format_moment_eqs` makes it visually more appealing: simplifies the symbolic
expressions further (expanding all terms), makes the time-dependence of moment
variables implicit (removes all `(t)`) and removes all trailing zeros (`2.0` → `2`).
"""
function format_moment_eqs(eqs::MomentEquations)

    odes = eqs.odes
    exprs  = []
    for i in 1:size(odes.eqs)[1]
        key = odes.states[i]
        eq = odes.eqs[i].rhs
        eq = isa(eq, Number) ? eq : clean_expr(eq)
        expr = "d"*string(key)*"/dt = "*string(eq)
        expr = replace(expr, "(t)"=>"")
        expr = replace(expr, ".0"=>"")
        push!(exprs, expr)
    end
    exprs

end

"""
    format_closure(eqs::ClosedMomentEquations)

Given [`ClosedMomentEquations`](@ref), return an array of formatted strings
representing closure functions for each higher order moment. The symbolic
expressions are formatted similarly to [`format_moment_eqs`](@ref) and can
be visualised using [Latexify](https://github.com/korsbo/Latexify.jl).
"""
function format_closure(eqs::ClosedMomentEquations)
    closure = eqs.closure
    exprs = []
    for i in keys(closure)
        eq = closure[i]
        eq = isa(eq, Number) ? eq : clean_expr(eq)
        expr = string(i)*" = "*string(eq)
        expr = replace(expr, "(t)"=>"")
        expr = replace(expr, ".0"=>"")
        push!(exprs, expr)
    end
    exprs
end


@latexrecipe function f(eqs::MomentEquations, type=:equations; inds::Array)

    env --> :align
    starred --> true
    cdot --> false

    if type == :equations
        return format_moment_eqs(eqs)
    elseif type == :closure
        return format_closure(eqs)
    else
        error("supported arguments are only `:equations` or `:closure`")
    end

end
