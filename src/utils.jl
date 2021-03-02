function sample_raw_moments(sol::EnsembleSolution, order::Int; naive::Bool=true)
    # need to build a t x n matrix for each timepoint, where t is the number of realisations (trajectories)
    # and n is the number of  variables

    # naive algorithm (naive=true), naivemoment(), performs much better in case of small N
    # when N is large, big performance improvement is expected(?) using moment() (naive=false)
    # this implementatation is experimental, need more testing to properly assert any
    # statement about performance + understand why the naive algorithm is still very fast
    # compared to the most trivial algorithm

    # TODO: rewrite this as sample_central_moments right now appearst to be faster

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


function deterministic_IC(μ₀::Vector, sys::MomentEquations, closed_sys::ODESystem)

    # sys as an argument is not strictly needed but it helps to implement checks that all arguments are consistent
    # closed_sys also needs to be included in case bernoulli variables were eliminated (so sys and closed_sys have a different no. of states)

    N = sys.N
    if N != length(μ₀)
        error("length of μ₀ and number of species are inconsistent")
    end
    no_states = length(closed_sys.states)

    μ_map = [sys.μ[iter] => μ₀[i] for (i, iter) in enumerate(sys.iter_1)]
    if typeof(sys) == CentralMomentEquations
        if string(closed_sys.states[N+1])[1] != 'M'
            error("the central moment equations and their closed counterparts are inconsistent")
        end
        moment_map = [closed_sys.states[i] => 0.0 for i in N+1:no_states]
    else
        if string(closed_sys.states[N+1])[1] != 'μ'
            error("the raw moment equations and their closed counterparts are inconsistent")
        end
        reverse_μ = Dict(iter => μ for (μ, iter) in sys.μ)
        #moment_map = []
        #for i in sys.N+1:length(closed_sys.states)
        #    iter = reverse_μ[closed_sys.states[i]]
        #    push!(moment_map, closed_sys.states[i] => prod(μ₀ .^ iter))
        #end
        moment_map = [closed_sys.states[i] => prod(μ₀ .^ reverse_μ[closed_sys.states[i]]) for i in N+1:no_states]
    end

    return vcat(μ_map, moment_map)

end

"""
    format_moment_eqs(odes::ODESystem)

Given a [`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/)
containing the moment equations, return an array of formatted strings representing
the ODEs that can be visualised using [Latexify](https://github.com/korsbo/Latexify.jl).
Although Latexify can be applied directly on the `ODESystem`, `format_moment_eqs`
makes it visually more appealing: simplifies the symbolic expressions further
(expanding all terms), makes the time-dependence of all moment variables implicit
(removes all `(t)`) and removes all trailing zeros (`2.0` → `2`).
"""
function format_moment_eqs(odes::ODESystem)

    # Input argument is ODESystem
    # Latexify applied on ODESystem directly introduces irrelevant multipliers (bug?)
    # for example, μ₁ may be turned into 1μ₁...
    # Here we reformat each ODE into a string expression and also remove the explicit
    # time dependence "(t)" of each variable as well as floats' trailing zeros".0"
    # then Latexify can be applied directly:
    # latexify(exprs, cdot=false, env=:align/:mdtable)
    # or latexify(exprs[1], cdot=false, env=:equation)
    # or the expressions can be copypasted
    # into Mathematica for example for further manipulation
    exprs  = []
    for i in 1:size(odes.eqs)[1]
        key = odes.states[i]
        eq = clean_expr(odes.eqs[i].rhs)
        expr = "d"*string(key)*"/dt = "*string(eq)
        expr = replace(expr, "(t)"=>"")
        expr = replace(expr, ".0"=>"")
        push!(exprs, expr)
    end
    exprs

end


function format_closure(closure::Dict{Any,Any})
    exprs = []
    for i in keys(closure)
        eq = clean_expr(closure[i])
        expr = string(i)*" = "*string(eq)
        expr = replace(expr, "(t)"=>"")
        expr = replace(expr, ".0"=>"")
        push!(exprs, expr)
    end
    exprs
end
