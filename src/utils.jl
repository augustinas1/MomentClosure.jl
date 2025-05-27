"""
    get_raw_moments(sol::EnsembleSolution, order::Int; naive::Bool=true, b::Int=2)

Given an `EnsembleSolution` of [DifferentialEquations ensemble simulation]
(https://diffeq.sciml.ai/stable/features/ensemble/#Performing-an-Ensemble-Simulation),
return a Dictionary of raw moments computed up to the specified `order` at each time step.

# Notes
- For example, the dictionary key `(2,0,1)` maps to an array containing the values of
  the raw moment ``μ_{201}`` at each time step.
- It is assumed that the time steps are all at the same time point for all trajectories
  (i.e., fixed `dt` used by the integrator or values were saved using `saveat`, as discussed [here]
  (https://diffeq.sciml.ai/stable/features/ensemble/#Time-steps-vs-time-points)).
- Moments are computed using [Cumulants.jl](https://github.com/iitis/Cumulants.jl) internally.
  The naive algorithm of moment tensor calculations (`naive=true`) is usually faster for small
  systems but the [proposed novel algorithm](https://arxiv.org/pdf/1701.05420.pdf) (`naive=false`)
  *should* be more efficient in case of many marginal variables. The [block size]
  (https://github.com/iitis/Cumulants.jl#block-size) `b` can also be specified for the
  novel algorithm and may have a significant effect on its performance.
- Only useful if higher order moments are needed: DifferentialEquations
  has a number of far more efficient and flexible ensemble statistics functions for
  means, variances and correlations, see [this tutorial]
  (https://diffeq.sciml.ai/stable/features/ensemble/#Analyzing-an-Ensemble-Experiment)
  for more details.
"""
function get_raw_moments(sol::EnsembleSolution, order::Int; naive::Bool=true, b::Int=2)

    # TODO: improve performance (it's rather poor now)
    # TODO: consider dropping the Dict and using DiffEqArrays

    if naive
        f_alg = naivemoment
        b = ()
    else
        f_alg = moment
    end

    no_t_pts = length(sol.u[1])
    N = length(sol.u[1].u[1])

    iter_all = construct_iter_all(N, order)
    iter_order = [filter(x -> sum(x) == m, iter_all) for m in 1:order]
    iter_μ = vcat(iter_order...)

    keys = Dict([iter => vcat([fill(i, n) for (i,n) in enumerate(iter)]...) for iter in iter_μ])
    μ = Dict([iter => Array{Float64}(undef, no_t_pts) for iter in iter_μ])

    for t_pt in 1:no_t_pts

        tslice = componentwise_vectors_timestep(sol, t_pt)
        tslice = hcat(tslice...) / 1

        for m in 1:order
            data = Array(f_alg(tslice, m, b...))
            for iter in iter_order[m]
                μ[iter][t_pt] = data[keys[iter]...]
            end
        end

    end

    μ

end

"""
    get_central_moments(sol::EnsembleSolution, order::Int; naive::Bool=true, b::Int=2)

Given an `EnsembleSolution` of [DifferentialEquations ensemble simulation]
(https://diffeq.sciml.ai/stable/features/ensemble/#Performing-an-Ensemble-Simulation),
return a Dictionary of central moment estimates computed up to the specified `order`
at each time step. See the notes of [`get_raw_moments`](@ref) function for more information.
"""
function get_central_moments(sol::EnsembleSolution, order::Int; naive::Bool=true, b::Int=2)

    if naive
        f_alg = naivemoment
        b = ()
    else
        f_alg = moment
    end

    no_t_pts = length(sol.u[1])
    N = length(sol.u[1].u[1])

    iter_all = construct_iter_all(N, order)
    iter_order = [filter(x-> sum(x) == m, iter_all) for m in 1:order]
    iter_M = filter(x -> sum(x) > 1, iter_all)

    keys = Dict([iter => vcat([fill(i, n) for (i,n) in enumerate(iter)]...) for iter in iter_all])
    M = Dict([iter => Array{Float64}(undef, no_t_pts) for iter in iter_M])

    μ = Dict()
    μ[Tuple(fill(0, N))] = 1.

    for t_pt in 1:no_t_pts
        tslice = componentwise_vectors_timestep(sol, t_pt)
        tslice = hcat(tslice...) / 1

        for m in 1:order
            data = Array(f_alg(tslice, m, b...))
            for iter in iter_order[m]
                μ[iter] = data[keys[iter]...]
            end
        end
        M_temp = raw_to_central_moments(N, order; μ)
        for iter in iter_M
            M[iter][t_pt] = M_temp[iter]
        end
    end

    M

end

"""
    get_cumulants(sol::EnsembleSolution, order::Int; naive::Bool=true, b::Int=2)

Given an `EnsembleSolution` of [DifferentialEquations ensemble simulation]
(https://diffeq.sciml.ai/stable/features/ensemble/#Performing-an-Ensemble-Simulation),
return a Dictionary of cumulant estimates computed up to the specified `order`
at each time step. See the notes of [`get_raw_moments`](@ref) function for more information.
"""
function get_cumulants(sol::EnsembleSolution, order::Int; naive::Bool=true, b::Int=2)

    no_t_pts = length(sol.u[1])
    N = length(sol.u[1].u[1])

    iter_all = construct_iter_all(N, order)
    iter_order = [filter(x-> sum(x) == m, iter_all) for m in 1:order]
    iter_κ = vcat(iter_order...)

    keys = Dict([iter => vcat([fill(i, n) for (i,n) in enumerate(iter)]...) for iter in iter_κ])
    κ = Dict([iter => Array{Float64}(undef, no_t_pts) for iter in iter_κ])

    if naive

        for t_pt in 1:no_t_pts

            tslice = componentwise_vectors_timestep(sol, t_pt)
            tslice = hcat(tslice...) / 1

            for m in 1:order
                data = Array(naivecumulant(tslice, m))
                for iter in iter_order[m]
                    κ[iter][t_pt] = data[keys[iter]...]
                end
            end

        end

    else

        for t_pt in 1:no_t_pts

            tslice = componentwise_vectors_timestep(sol, t_pt)
            tslice = hcat(tslice...) / 1
            κ_tensor = cumulants(tslice, order, b)

            for m in 1:order

                data = Array(κ_tensor[m])

                for iter in iter_order[m]
                    κ[iter][t_pt] = data[keys[iter]...]
                end

            end

        end

    end

    κ

end


"""
    get_moments_FSP(sol::ODESolution, order::Int, moment_type::String)

Given an `ODESolution` obtained using [FiniteStateProjection.jl](https://github.com/kaandocal/FiniteStateProjection.jl),
return a Dictionary of moments computed up to the specified `order` at each time step.
Here, `moment_type` specifies the type of moments to be computed: available options are `raw`,
`central` or `cumulant`.

# Notes
- The `ODESolution` represents the time-evolution of the probability density function
  that is the solution of the Chemical Master Equation approximated using Finite State Projection
  algorithms. See the documentation of [FiniteStateProjection.jl](https://github.com/kaandocal/FiniteStateProjection.jl)
  for more information.
"""
function get_moments_FSP(sol::ODESolution, order::Int, moment_type::String)

    @assert any(moment_type .== ["raw", "central", "cumulant"]) "Moment type "*moment_type*" does not exist"

    state_space = size(sol)[1:end-1]
    N = length(state_space)
    no_t_pts = length(sol.u)

    iter_moments = MomentClosure.construct_iter_all(N, order)[2:end]
    moments = Dict([iter => Array{Float64}(undef, no_t_pts) for iter in iter_moments])

    iter_state = Iterators.product((0:i-1 for i in state_space)...)

    if moment_type == "raw"
        for t_pt in 1:no_t_pts
            tslice = sol[t_pt]
            for μ_ind in iter_moments
                moments[μ_ind][t_pt] = sum(prod(n_state.^μ_ind) * tslice[(n_state.+1)...] for n_state in iter_state)
            end
        end
    else
        μ = Dict{NTuple{N, Float64}}()
        μ[Tuple(zeros(N))] = 1.

        for t_pt in 1:no_t_pts
            tslice = sol[t_pt]
            for μ_ind in iter_moments
                μ[μ_ind] = sum(prod(n_state.^μ_ind) * tslice[(n_state.+1)...] for n_state in iter_state)
            end

            if moment_type == "central"
                moment_temp = MomentClosure.raw_to_central_moments(N, order; μ)
            else
                moment_temp = MomentClosure.cumulants_to_raw_moments(N, order; μ)
            end

            for iter in iter_moments
                moments[iter][t_pt] = moment_temp[iter]
            end
        end
    end

    moments

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
  in the corresponding reaction system (can be checked with the 
  [`Catalyst.speciesmap`](https://docs.sciml.ai/Catalyst/stable/api/core_api/#Catalyst.speciesmap) 
  function).
- As higher-order moment functions under log-normal, gamma, derivative matching and
  the conditional closures involve moments raised to negative powers, setting initial
  molecule numbers of certain species to *zeros* will result in NaN errors when solving
  the ODEs (the specifics depend on the system at hand).
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

    vars = unknowns(odes)
    no_states = length(vars)
    if typeof(sys) == CentralMomentEquations
        moment_map = [vars[i] => 0.0 for i in N+1:no_states]
    else
        reverse_μ = Dict(μ => iter for (iter, μ) in sys.μ)
        moment_map = [vars[i] => prod(u₀ .^ reverse_μ[vars[i]]) for i in N+1:no_states]
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

    sys = eqs.odes
    odes = get_eqs(sys)
    vars = unknowns(sys)
    exprs  = []

    for i in 1:size(odes)[1]
        key = vars[i]
        eq = odes[i].rhs
        eq = expand_mod(eq)
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

# Notes
- `include_all` argument specifies whether all higher-order moments (`include_all=true`)
  or only those encountered in the moment ODEs specifically (`include_all=false`,
  the default option) will be returned.
"""
function format_closure(eqs::ClosedMomentEquations; format_all::Bool=false)
    closure = eqs.closure
    exprs = []

    if format_all
        iter = keys(closure)
    else
        iter = setdiff(unknowns(eqs.open_eqs.odes), unknowns(eqs.odes))
    end

    for i in iter
        eq = closure[i]
        eq = expand_mod(eq)
        expr = string(i)*" = "*string(eq)
        expr = replace(expr, "(t)"=>"")
        expr = replace(expr, ".0"=>"")
        push!(exprs, expr)
    end
    exprs
end


@latexrecipe function f(eqs::MomentEquations, type=:equations; print_all=false)

    env --> :align
    starred --> true
    cdot --> false

    if type == :equations
        return format_moment_eqs(eqs)
    elseif type == :closure
        return format_closure(eqs, format_all=print_all)
    else
        error("supported arguments are only `:equations` or `:closure`")
    end

end
