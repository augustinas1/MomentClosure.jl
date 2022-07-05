using MomentClosure
using MomentClosure: gen_iter
using SciMLBase.EnsembleAnalysis, JumpProcesses
using Test

# default SIR model example from JumpProcesses
substoich = [[1 => 1, 2 => 1],
             [2 => 1]]        
netstoich = [[1 => -1, 2 => 1],
             [2 => -1, 3 => 1]]
p     = (0.1/1000, 0.01)
pidxs = [1, 2]
maj   = MassActionJump(substoich, netstoich; param_idxs=pidxs)
u₀    = [9, 1, 0] 
tspan = (0.0, 50.0)
dprob = DiscreteProblem(u₀, tspan, p)
jprob = JumpProblem(dprob, Direct(), maj, save_positions=(false, false))
@time sol   = solve(EnsembleProblem(jprob), saveat=25., SSAStepper(), trajectories=100)

μ_1 = get_raw_moments(sol, 2)
μ_2 = get_raw_moments(sol, 2, naive=false, b=3)
@test isequal(μ_1, μ_2)

M_1 = get_central_moments(sol, 2)
M_2 = get_central_moments(sol, 2, naive=false, b=3)
@test isequal(M_1, M_2)

κ_1 = get_cumulants(sol, 2)
κ_2 = get_cumulants(sol, 2, naive=false, b=3)
@test isequal(κ_1, κ_2)

m_SSA, var_SSA = timeseries_steps_meanvar(sol)
μ = μ_1; M = M_1; κ = κ_1
inds = gen_iter(3, 1)
@test all(isapprox(μ[inds[i]], m_SSA[i, :]) for i in 1:3)

inds = [ind .* 2 for ind in inds]
@test all(isapprox(M[inds[i]], var_SSA[i, :], atol=1e-2) for i in 1:3)

@test isapprox(μ[(1, 1, 0)] .- μ[(0, 1, 0)].*μ[(1,0,0)], M[(1, 1, 0)])
@test isapprox(μ[(1, 0, 1)] .- μ[(1, 0, 0)].*μ[(0,0,1)], M[(1, 0, 1)])

inds = gen_iter(3, 1)
@test all(isapprox(μ[i], κ[i]) for i in gen_iter(3,1))
@test all(isapprox(M[i], κ[i]) for i in gen_iter(3, 2))