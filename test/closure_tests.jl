using MomentClosure
using MomentClosure: value
using Test, SafeTestsets
using Catalyst

@testset "moment closure tests" begin

    @parameters t, c₁, c₂, c₃, c₄, Ω
    @variables X(t), Y(t)

    S_mat = [ 1 -1  1 -1;
             -1  1  0  0]
    a = [c₁*X*Y*(X-1)/Ω^2, c₂*X, c₃*Ω, c₄*X]
    rn = ReactionSystemMod(t, [X, Y], [c₁, c₂, c₃, c₄, Ω], a, S_mat)

    sys = generate_central_moment_eqs(rn, 2, 4)
    expr1 = sys.odes.eqs[1].rhs
    closed_odes, closure, closure_symbolic = moment_closure(sys, "zero")
    @test closure[sys.M[2,2]] == 0
    expr1 = closed_odes.eqs[1].rhs
    M = sys.M
    μ = sys.μ
    expr2 = c₃*Ω + M[1,1]*c₁*μ[1,0]*(Ω^-2) + M[1,1]*c₁*(Ω^-2)*(μ[1,0]- 1) + c₁*M[2,0]*μ[0,1]*(Ω^-2) +
            c₁*μ[0,1]*μ[1,0]*(Ω^-2)*(μ[1,0] - 1) - c₂*μ[1,0] - c₄*μ[1,0]
    expr2 = simplify(value.(expr2))
    @test isequal(expand_mod(expr1), expand_mod(expr2))

    sys = generate_raw_moment_eqs(rn, 2)
    closed_odes, closure, closure_symbolic = moment_closure(sys, "zero")
    μ = sys.μ
    expr1 = closure[μ[3,0]]
    expr2 = -2*μ[1,0]^3 + 3*μ[1,0]*μ[2,0]
    @test isequal(expr1, expr2)

end
