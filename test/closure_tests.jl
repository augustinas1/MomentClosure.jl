using MomentClosure
using MomentClosure: value, define_M, define_μ
using Test
using Catalyst

@parameters t, c₁, c₂, c₃, c₄, Ω
@variables X(t), Y(t)

S_mat = [ 1 -1  1 -1;
         -1  1  0  0]
a = [c₁*X*Y*(X-1)/Ω^2, c₂*X, c₃*Ω, c₄*X]
rn = ReactionSystemMod(t, [X, Y], [c₁, c₂, c₃, c₄, Ω], a, S_mat)

μ = define_μ(2,4)
M = define_M(2,4)
sys = generate_central_moment_eqs(rn, 2, 4)
expr1 = sys.odes.eqs[1].rhs

closed_odes, closure = moment_closure(sys, "zero")
@test closure[M[2,2]] == 0
expr1 = closed_odes.eqs[1].rhs
expr2 = c₃*Ω + M[1,1]*c₁*μ[1,0]*(Ω^-2) + M[1,1]*c₁*(Ω^-2)*(μ[1,0]- 1) + c₁*M[2,0]*μ[0,1]*(Ω^-2) +
        c₁*μ[0,1]*μ[1,0]*(Ω^-2)*(μ[1,0] - 1) - c₂*μ[1,0] - c₄*μ[1,0]
expr2 = simplify(value.(expr2))
@test isequal(expand_mod(expr1), expand_mod(expr2))

closed_odes, closure = moment_closure(sys, "normal")
@test isequal(closure[M[0,4]], 3*M[0,2]^2)

closed_odes, closure = moment_closure(sys, "log-normal")
expr1 = closure[M[1,2]]
expr2 = exp(log(μ[1,0]) + log(1 + M[0,2]*μ[0,1]^-2) + 2(log(μ[0,1]) + log(1 + M[1,1]*(μ[0,1]^-1)*(μ[1,0]^-1)))) -
    M[0,2]*μ[1,0] - μ[1,0]*μ[0,1]^2 - 2*M[1,1]*μ[0,1]
@test isequal(expr1, expr2)

closed_odes, closure = moment_closure(sys, "poisson")
@test isequal(closure[M[3,0]], μ[1,0])

closed_odes, closure = moment_closure(sys, "gamma")
expr1 = closure[M[2,1]]
expr2 = M[2,0]*μ[0,1] + μ[0,1]*μ[1,0]^2 + 2*M[1,1]*μ[1,0] + 2*M[1,1]*M[2,0]*μ[1,0]^-1 -
    M[2,0]*μ[0,1] - μ[0,1]*μ[1,0]^2 - 2*M[1,1]*μ[1,0]
@test isequal(expr1, expr2)

closed_odes, closure = moment_closure(sys, "derivative matching")
expr1 = closure[sys.M[0,4]]
expr2 = μ[0,1]^4*(M[0,2]+μ[0,1]^2)^-6*(M[0,3]+μ[0,1]^3+3*M[0,2]*μ[0,1])^4 - μ[0,1]^4 -
    6*M[0,2]*μ[0,1]^2 - 4*M[0,3]*μ[0,1]
@test isequal(expr1, expr2)

sys = generate_raw_moment_eqs(rn, 2)
closed_odes, closure = moment_closure(sys, "zero")
expr1 = closure[μ[3,0]]
expr2 = -2*μ[1,0]^3 + 3*μ[1,0]*μ[2,0]
@test isequal(expr1, expr2)

closed_odes, closure = moment_closure(sys, "normal")
expr1 = closure[sys.μ[1,2]]
expr2 = 2*μ[0,1]*μ[1,1] + μ[0,2]*μ[1,0] - 2*μ[1,0]*μ[0,1]^2
@test isequal(expr1, expr2)
expr1 = closed_odes.eqs[5].rhs
expr2 = c₂*μ[1,0] + 2*c₂*μ[1,1] + c₁*μ[0,1]*μ[2,0]/Ω^2 - c₁*μ[1,1]/Ω^2 - 4*c₁*μ[1,1]^2/Ω^2 + 4*c₁*μ[0,1]*μ[1,1]/Ω^2 -
    2*c₁*μ[0,1]*μ[1,0]^2/Ω^2 + 2*c₁*μ[0,2]*μ[1,0]/Ω^2 -2*c₁*μ[0,2]*μ[2,0]/Ω^2 + 2*c₁*μ[1,0]*μ[1,1]/Ω^2 -
    4*c₁*μ[1,0]*μ[0,1]^2/Ω^2 + 4*c₁*μ[0,1]^2*μ[1,0]^2/Ω^2
expr2 = simplify(value.(expr2))
@test isequal(expand_mod(expr1), expr2)

closed_odes, closure = moment_closure(sys, "log-normal")
expr1 = closure[μ[1,3]]
expr2 = exp(log(μ[1,0]) + 3*(log(μ[0,2]*μ[0,1]^(-2)) + log(μ[1,1]*(μ[0,1]^-1)*(μ[1,0]^-1)) + log(μ[0,1])))
@test isequal(expr1, expr2)

closed_odes, closure = moment_closure(sys, "poisson")
expr1 = closure[μ[4,0]]
expr2 = μ[1,0] + 6 * μ[1,0]^4 + 3*μ[2,0]^2 + 4*μ[1,0]*μ[3,0] -12*μ[2,0]*μ[1,0]^2
@test isequal(expr1, expr2)

closed_odes, closure = moment_closure(sys, "gamma")
expr1 = closure[μ[0,3]]
expr2 = μ[0,1]^3 + 3*μ[0,1]*(μ[0,2]-μ[0,1]^2) + 2*μ[0,1]^-1*(μ[0,2]-μ[0,1]^2)^2
@test isequal(expr1, expr2)

closed_odes, closure = moment_closure(sys, "derivative matching")
expr1 = closure[sys.μ[0,3]]
expr2 = μ[0,1]^-3*μ[0,2]^3
@test isequal(expr1, expr2)
