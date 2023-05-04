using MomentClosure
using MomentClosure: define_M, define_μ
using Symbolics: value, expand
using Test
using Catalyst

rn = @reaction_network begin
    (c₁/Ω^2), 2X + Y → 3X
    (c₂), X → Y
    (Ω*c₃, c₄), 0 ↔ X
end

@syms c₁::Real c₂::Real c₃::Real c₄::Real Ω::Real
μ = define_μ(2,4)
M = define_M(2,4)

# --- Test closures on central moment equations ---

sys = generate_central_moment_eqs(rn, 2, 4, combinatoric_ratelaws=false)
expr1 = sys.odes.eqs[1].rhs
closed_eqs = moment_closure(sys, "zero")
@test closed_eqs.closure[M[2,2]] == 0
expr1 = closed_eqs.odes.eqs[1].rhs
expr2 = c₃*Ω + M[1,1]*c₁*μ[1,0]*(Ω^-2) + M[1,1]*c₁*(Ω^-2)*(μ[1,0]- 1) + c₁*M[2,0]*μ[0,1]*(Ω^-2) +
        c₁*μ[0,1]*μ[1,0]*(Ω^-2)*(μ[1,0] - 1) - c₂*μ[1,0] - c₄*μ[1,0]
@test isequal(simplify(expr1), simplify(expr2))

# check that deterministic_IC is working with central moments
ic_values = Dict(deterministic_IC([2, 5], closed_eqs))
@test length(ic_values) == 5
@test ic_values[μ[1,0]] == 2 && ic_values[μ[0,1]] == 5 && ic_values[M[1,1]] == 0

closed_eqs = moment_closure(sys, "normal")
@test isequal(closed_eqs.closure[M[0,4]], 3*M[0,2]^2)

closed_eqs= moment_closure(sys, "log-normal")
expr1 = closed_eqs.closure[M[1,2]]
expr2 = μ[1,0]*μ[0,1]^2*(1.0+M[0,2]*μ[0,1]^-2)*(1.0 + M[1,1]*(μ[0,1]^-1)*(μ[1,0]^-1))^2 -
        M[0,2]*μ[1,0] - μ[1,0]*μ[0,1]^2 - 2*M[1,1]*μ[0,1]
@test isequal(expr1, simplify(expand(expr2)))

closed_eqs = moment_closure(sys, "poisson")
@test isequal(closed_eqs.closure[M[3,0]], μ[1,0])

closed_eqs = moment_closure(sys, "gamma")
expr1 = closed_eqs.closure[M[2,1]]
expr2 = M[2,0]*μ[0,1] + μ[0,1]*μ[1,0]^2 + 2*M[1,1]*μ[1,0] + 2*M[1,1]*M[2,0]*μ[1,0]^-1 -
    M[2,0]*μ[0,1] - μ[0,1]*μ[1,0]^2 - 2*M[1,1]*μ[1,0]
@test isequal(expr1, expr2)

closed_eqs = moment_closure(sys, "derivative matching")
expr1 = closed_eqs.closure[sys.M[0,4]]
expr2 = μ[0,1]^4*(M[0,2]+μ[0,1]^2)^-6*(M[0,3]+μ[0,1]^3+3*M[0,2]*μ[0,1])^4 - μ[0,1]^4 -
    6*M[0,2]*μ[0,1]^2 - 4*M[0,3]*μ[0,1]
@test isequal(expand(expr1), expand(simplify(expr2)))

# --- Test closures on raw moment equations ---

sys = generate_raw_moment_eqs(rn, 2, combinatoric_ratelaws=false)
closed_eqs = moment_closure(sys, "zero")
expr1 = closed_eqs.closure[μ[3,0]]
expr2 = -2*μ[1,0]^3 + 3*μ[1,0]*μ[2,0]
@test isequal(expr1, expr2)

closed_eqs = moment_closure(sys, "normal")
expr1 = closed_eqs.closure[sys.μ[1,2]]
expr2 = 2*μ[0,1]*μ[1,1] + μ[0,2]*μ[1,0] - 2*μ[1,0]*μ[0,1]^2
@test isequal(expr1, expr2)
expr1 = closed_eqs.odes.eqs[5].rhs
expr2 = c₂*μ[1,0] + 2*c₂*μ[1,1] + c₁*μ[0,1]*μ[2,0]*Ω^-2 - c₁*μ[1,1]*Ω^-2 - 4*c₁*μ[1,1]^2*Ω^-2 + 4*c₁*μ[0,1]*μ[1,1]*Ω^-2 -
    2*c₁*μ[0,1]*μ[1,0]^2*Ω^-2 + 2*c₁*μ[0,2]*μ[1,0]*Ω^-2 -2*c₁*μ[0,2]*μ[2,0]*Ω^-2 + 2*c₁*μ[1,0]*μ[1,1]*Ω^-2 -
    4*c₁*μ[1,0]*μ[0,1]^2*Ω^-2 + 4*c₁*μ[0,1]^2*μ[1,0]^2*Ω^-2
@test isequal(simplify(expr1), simplify(expr2))

# check that deterministic_IC is working with raw moments
ic_values = Dict(deterministic_IC([2, 5], closed_eqs))
@test length(ic_values) == 5
@test ic_values[μ[1,0]] == 2 && ic_values[μ[0,1]] == 5 && ic_values[μ[1,1]] == 10

closed_eqs = moment_closure(sys, "log-normal")
expr1 = closed_eqs.closure[μ[1,3]]
expr2 = μ[0,1]^-6 * μ[0,2]^3 * μ[1,0]^-2 * μ[1,1]^3
@test isequal(expr1, expr2)

closed_eqs = moment_closure(sys, "poisson")
expr1 = closed_eqs.closure[μ[4,0]]
expr2 = μ[1,0] + 6 * μ[1,0]^4 + 3*μ[2,0]^2 + 4*μ[1,0]*μ[3,0] -12*μ[2,0]*μ[1,0]^2
@test isequal(expr1, expr2)

closed_eqs = moment_closure(sys, "gamma")
expr1 = closed_eqs.closure[μ[0,3]]
expr2 = μ[0,1]^3 + 3*μ[0,1]*(μ[0,2]-μ[0,1]^2) + 2*μ[0,1]^-1*(μ[0,2]-μ[0,1]^2)^2
@test isequal(simplify(expr1), simplify(expr2))

closed_eqs = moment_closure(sys, "derivative matching")
expr1 = closed_eqs.closure[sys.μ[0,3]]
expr2 = μ[0,1]^-3*μ[0,2]^3
@test isequal(expr1, expr2)
