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
    (k₃ * x * y / (x + k₇)), x ⇒ 0
    (k₄ * x), 0 ⇒ y₀
    (k₅), y₀ → y
    (k₆), y → 0
end

eqs = generate_central_moment_eqs(rn, 2, 3, combinatoric_ratelaws = false)
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
# Reference is the converged solution (stable to ~10 digits for reltol/abstol ≤ 1e-10);
# integrate at a tight tolerance so the comparison is independent of solver defaults.
sol = solve(oprob, Tsit5(), saveat = 10:10:50, reltol = 1.0e-10, abstol = 1.0e-10)
res = [35.156249675, 60.255177096, 29.237221341, 65.59590462, 24.900095364]
@test isapprox(sol[1, :], res, atol = 1.0e-6)
oprob = ODEProblem(closed_eqs, oprob.u0, tspan, pmap, use_deterministic_IC = false)
sol = solve(oprob, Tsit5(), saveat = 10:10:50, reltol = 1.0e-10, abstol = 1.0e-10)
@test isapprox(sol[1, :], res, atol = 1.0e-6)
