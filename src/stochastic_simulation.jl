#=
    JumpProblem(rn::ReactionSystemMod, prob::DiscreteProblem, aggregator, args...; kwargs...)

Given a `ReactionSystemMod`, DifferentialEquations [`DiscreteProblem`]
(https://diffeq.sciml.ai/stable/types/discrete_types/#Discrete-Problems) and [the method
for aggregating constant jumps](https://diffeq.sciml.ai/stable/types/jump_types/#Jump-Problems),
return a corresponding  [`JumpProblem`](https://diffeq.sciml.ai/stable/types/jump_types/#Jump-Problems)
which can be solved using standard DifferentialEquations methods.

# Notes
- This function essentially massages a `ReactionSystemMod` into a form which could be used
  with DifferentialEquations SSA solvers. Namely, a `ReactionSystemMod` is converted into a
  system of [ConstantRateJumps]
  (https://diffeq.sciml.ai/stable/types/jump_types/#Types-of-Jumps:-Regular,-Variable,-Constant-Rate-and-Mass-Action)
  taking into account the geometrically distributed reaction products.
- The implementation here is very restricted: time-dependent propensities are currently not
  allowed; [MassActionJumps]
  (https://diffeq.sciml.ai/stable/types/jump_types/#Types-of-Jumps:-Regular,-Variable,-Constant-Rate-and-Mass-Action)
  cannot be used due to limitations of the `ReactionSystemMod` API, hence the SSA performance is subpar.
=#
function DiffEqJump.JumpProblem(rn::ReactionSystemMod, prob::DiscreteProblem, aggregator, args...; kwargs...)

    # TODO:
    # - add random seed for sampling from geometric
    # - allow VariableRateJumps
    # - using MassActionJumps where possible would be significantly faster but they require
    #  the reactant stoichometry matrix and knowing which reactions are mass-action (not possible with current API)

    # construct the mapping from each paramater to their value if the DiscreteProblem
    # contains a vector of parameter values
    if prob.p isa Array{<:Real, 1}
        pmap = Dict(Pair.(reactionparams(rn), prob.p / 1.))
    else
        # otherwise we assume `prob` already contains the required mapping
        pmap = Dict(prob.p)
    end

    # build rate functions
    as = [substitute(a, pmap) for a in propensities(rn)]
    # in case of VariableRateJumps or symbolic params this should be amended
    af = [eval(build_function(a, species(rn)...)) for a in as]

    S_mat = netstoichmat(rn)

    is_geometric = fill(false, size(S_mat))
    rngs = similar(S_mat, Any)

    for (i, j) in Tuple.(findall(x->isa(x, Sym), S_mat))
        is_geometric[i, j] = true
        m = pmap[S_mat[i, j]]
        # note that Geometric considers success probability p
        # the mean is (1-p)/p, hence p = 1/(mean+1)
        rngs[i, j] = Geometric(1/(m+1))
    end

    jumps = []

    for r in 1:numreactions(rn)

        rate(u, p, t) = af[r](u...)

        if sum(is_geometric[:, r]) > 0

            inds_n = findall(s->!s, is_geometric[:, r])
            inds_g = findall(is_geometric[:, r])
            function affect1!(integrator)
                integrator.u[inds_n] .+= S_mat[inds_n, r]
                for ind in inds_g
                    integrator.u[ind] += rand(rngs[ind, r])
                end
            end
            push!(jumps, ConstantRateJump(rate, affect1!))
        else

            function affect2!(integrator)
                integrator.u .+= S_mat[:, r]
            end
            push!(jumps, ConstantRateJump(rate, affect2!))
        end


    end

    return JumpProblem(prob, aggregator, jumps..., args...; kwargs...)

end
