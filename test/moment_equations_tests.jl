using MomentClosure
using MomentClosure: expected_coeff
using ModelingToolkit: value, simplify
using SymbolicUtils: expand
using Catalyst
using Test

@parameters b
b = b.val
@test expected_coeff(2, 3) == 8
@test isequal(expected_coeff(b, 3), b*(1 + (6b) + (6(b^2))))
expr = ((b*((1 + b)^-1)) + ((b^6)*((1 + b)^-6)) + (57(b^2)*((1 + b)^-2)) + (302(b^3)*((1 + b)^-3)) + (302(b^4)*((1 + b)^-4)) + (57(b^5)*((1 + b)^-5)))*((1 + b)^6)
@test isequal(expected_coeff(b, 6), expr)

@parameters t, c₁, c₂, c₃, c₄, Ω
@variables X(t), Y(t)

S_mat = [ 1 -1  1 -1;
         -1  1  0  0]
a = [c₁*X*Y*(X-1)/Ω^2, c₂*X, c₃*Ω, c₄*X]
rn = ReactionSystemMod(t, [X, Y], [c₁, c₂, c₃, c₄, Ω], a, S_mat)

sys = generate_central_moment_eqs(rn, 2, 4)
expr1 = sys.odes.eqs[2].rhs
@test isequal(MomentClosure.Differential(t)(sys.μ[1,0]), sys.odes.eqs[1].lhs)
μ = sys.μ
M = sys.M
# simplify(value()) is needed due to diffs between ModelingToolkit Num type and SymbolicUtils types
expr2 = c₂*μ[1,0] + c₁*M[1,1]*(Ω^-2) + c₁*μ[0,1]*μ[1,0]*Ω^-2 - c₁*M[2,1]*Ω^-2 -
        2*c₁*M[1,1]*μ[1,0]*Ω^-2 - c₁*M[2,0]*μ[0,1]*Ω^-2 - c₁*μ[0,1]*Ω^-2*μ[1,0]^2
expr2 = simplify(value.(expr2))
@test isequal(expand(expr1), expr2)

sys = generate_central_moment_eqs(rn, 2)
expr1 = sys.odes.eqs[2].rhs
@test isequal(MomentClosure.Differential(t)(sys.μ[1,0]), sys.odes.eqs[1].lhs)
@test isequal(expand(expr1), expr2)

sys = generate_raw_moment_eqs(rn, 2)
μ = sys.μ
expr1 = sys.odes.eqs[4].rhs
expr2 = c₂*μ[2,0] + c₁*μ[1,1]/Ω^2 - c₁*μ[3,1]/Ω^2 -c₂*(μ[1,0] + μ[1,1]) +
        c₃*μ[0,1]*Ω - c₄*μ[1,1] - c₁*μ[1,2]/Ω^2 + c₁*μ[2,2]/Ω^2
expr2 = simplify(value.(expr2))
@test isequal(expand(expr1), expand(expr2))

# corner case - a linear propensity with rate coefficient 1
rn = @reaction_network begin
    1, X → 0
end
sys = generate_raw_moment_eqs(rn, 2)
μ = sys.μ
@test isequal(MomentClosure.Differential(t)(μ[(1,)]) ~ -μ[(1,)], sys.odes.eqs[1])

# time-dependent propensity tests

rn = @reaction_network begin
        (X^1.5), X ⇒ Y
end c₁
@test_throws ErrorException generate_raw_moment_eqs(rn, 2)

rn = @reaction_network begin
        (c₁*Y^3*X^1.5), X ⇒ Y
end c₁
@test_throws ErrorException generate_raw_moment_eqs(rn, 2)

rn = @reaction_network begin
        (c₁*sin(exp(X)+Y)), X ⇒ Y
end c₁
@test_throws ErrorException generate_raw_moment_eqs(rn, 2)

rn = @reaction_network begin
        (c₁*sin(exp(t))+Y^(3+1)), X ⇒ Y
end c₁
eqs = generate_raw_moment_eqs(rn, 2)
@test eqs.q_order == 5
