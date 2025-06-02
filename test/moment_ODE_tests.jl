using MomentClosure
using Test
using Catalyst
using OrdinaryDiffEqTsit5
using ModelingToolkit: t_nounits as t

@parameters k₁ k₂ k₃ k₄ k₅ k₆ k₇
@species x(t) y(t) y₀(t)

rn = @reaction_network begin
    @parameters k₁ k₂ k₃ k₄ k₅ k₆ k₇
    (k₁), 0 → x
    (k₂), x → 0
    (k₃*x*y/(x+k₇)), x ⇒ 0
    (k₄*x), 0 ⇒ y₀
    (k₅), y₀ → y
    (k₆), y → 0
end

eqs = generate_central_moment_eqs(rn, 2, 3, combinatoric_ratelaws=false)
closed_eqs = moment_closure(eqs, "normal")

# test deterministic_IC
idxs = [speciesmap(rn)[x], speciesmap(rn)[y], speciesmap(rn)[y₀]]
u0map1 = (:y₀ => 30, :x => 70, :y => 60)
u0map2 = Dict(x => 70, y => 60, y₀ => 30)
u0map3 = [70, 60, 30][idxs]
@test isequal(deterministic_IC(u0map1, closed_eqs), deterministic_IC(u0map2, closed_eqs))
@test isequal(deterministic_IC(u0map1, closed_eqs), deterministic_IC(u0map3, closed_eqs))

# test ODEProblem
pmap = [:k₁ => 90, :k₂ => 0.002, :k₃ => 1.7, :k₄ => 1.1, :k₅ => 0.93, :k₆ => 0.96, :k₇ => 0.01]
tspan = (0, 200)
oprob = ODEProblem(closed_eqs, u0map1, tspan, pmap)
sol = solve(oprob, Tsit5(), saveat=10:10:50)
res = [35.155, 60.259, 29.232, 65.602, 24.891]
@test isapprox(sol[1, :], res, atol=1e-3)
oprob = ODEProblem(closed_eqs, oprob.u0, tspan, pmap, use_deterministic_IC=false)
sol = solve(oprob, Tsit5(), saveat=10:10:50)
@test isapprox(sol[1, :], res, atol=1e-3)