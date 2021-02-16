function compute_central_moment(sol, i_vec, mu=nothing)

    # trivial implementation of a biased multivariate central moment estimator
    # if calculating variances it is more efficient to use timeseries_steps_meanvar()
    # from DifferentialEquations.EnsembleAnalysis
    # NOTE: the function as currently implemented uses DifferentialEquations
    # Can be extremely slow...

    # TODO: is "EnsembleAnalysis." required if not explicitly loaded?

    inds = findall(x->x!=0, i_vec)
    i_vec = i_vec[inds]

    # NOTE: requires DifferentialEquations.EnsembleAnalysis
    if mu==nothing
        mu = EnsembleAnalysis.timeseries_steps_mean(sol)
    end

    M = zeros(length(sol[1]))
    for t in 1:length(sol[1])  # iterating over time steps
        for path in 1:length(sol) # iterating over trajectories
            M[t] += prod([(sol[path][t][idx] - mu[t][idx])^i for (idx, i) in zip(inds, i_vec)])
        end
        M[t] = M[t]/(length(sol))
    end

    return DiffEqArray(M, sol[1].t)

end


function compute_raw_moment(sol, i_vec, mu=nothing)

    # trivial implementation of a biased estimator of raw moments
    # if calculating means it is more efficient
    # to use timeseries_steps_meanvar() from DifferentialEquations.jl
    # TODO: remove mu

    inds = findall(x->x!=0, i_vec)
    i_vec = i_vec[inds]

    if mu==nothing
        mu = EnsembleAnalysis.timeseries_steps_mean(sol)
    end

    M = zeros(length(sol[1]))
    for t in 1:length(sol[1])  # iterating over time steps
        for path in 1:length(sol) # iterating over trajectories
            M[t] += prod([(sol[path][t][idx])^i for (idx, i) in zip(inds, i_vec)])
        end
        M[t] = M[t]/(length(sol))
    end

    return DiffEqArray(M, sol[1].t)

end


function deterministic_IC(μ₀::Vector, sys::Union{RawMomentEquations, CentralMomentEquations}, closed_sys::ODESystem)

    # sys as an argument is not strictly needed but it helps to implement checks that all arguments are consistent
    # closed_sys also needs to be included in case bernoulli variables were eliminated (so sys and closed_sys have a different no. of states)

    N = sys.N
    if N != length(μ₀)
        error("length of μ₀ and number of species are inconsistent")
    end
    no_states = length(closed_sys.states)

    μ_map = [sys.μ[iter] => μ₀[i] for (i, iter) in enumerate(sys.unit_vec)]
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
