using MomentClosure
using MomentClosure: expected_coeff, expand_div
using Symbolics: value, simplify, expand
using Distributions: Geometric
using Catalyst
using Test

# --- check that expected values of stoichiometric coefficients are computed correctly ---

@parameters b, d, p
b = value(b); d = value(d); p = value(p)
@test isequal(expected_coeff(2, 3), 8)
@test isequal(expected_coeff(b, 2), b^2)
@test isequal(expected_coeff(b*d, 2), b^2*d^2)
@test isequal(expected_coeff(b/d, 2), b^2/d^2)
@test isequal(expected_coeff(b+1, 2), b^2+2b+1)
@test isequal(expected_coeff(sin(b), 2), sin(b)^2)
@register_symbolic Geometric(p)
m = rand(Geometric(p))
@test isequal(expand_div(expected_coeff(m, 1)), p^-1 - 1)
@test isequal(expand_div(expected_coeff(m, 3)), expand(-1 + 7/p - 12/p^2 + 6/p^3))
@register_symbolic Bernoulli(p)
@test_throws AssertionError expected_coeff(rand(Bernoulli(p)), 1)

# --- check that moment equations are generated correctly ---

@parameters t, c₁, c₂, c₃, c₄, Ω

rn = @reaction_network begin
    (c₁*Ω^-2), 2X + Y → 3X
    (c₂), X → Y
    (Ω*c₃, c₄), 0 ↔ X
end

sys = generate_central_moment_eqs(rn, 2, 4, combinatoric_ratelaw=false)
expr1 = sys.odes.eqs[2].rhs
@test isequal(MomentClosure.Differential(t)(sys.μ[1,0]), sys.odes.eqs[1].lhs)
μ = sys.μ
M = sys.M
# simplify(value()) is needed due to diffs between ModelingToolkit Num type and SymbolicUtils types
expr2 = c₂*μ[1,0] + c₁*M[1,1]*(Ω^-2) + c₁*μ[0,1]*μ[1,0]*Ω^-2 - c₁*M[2,1]*Ω^-2 -
        2*c₁*M[1,1]*μ[1,0]*Ω^-2 - c₁*M[2,0]*μ[0,1]*Ω^-2 - c₁*μ[0,1]*Ω^-2*μ[1,0]^2
expr2 = simplify(value.(expr2))
#@test isequal(expand(expr1), expr2)
@test isequal(expr1, expr2)


sys = generate_central_moment_eqs(rn, 2, combinatoric_ratelaw=false)
expr1 = sys.odes.eqs[2].rhs
@test isequal(MomentClosure.Differential(t)(sys.μ[1,0]), sys.odes.eqs[1].lhs)
@test isequal(expr1, expr2)

sys = generate_raw_moment_eqs(rn, 2, combinatoric_ratelaw=false)
μ = sys.μ
expr1 = sys.odes.eqs[4].rhs
expr2 = c₂*μ[2,0] + c₁*μ[1,1]*Ω^-2 - c₁*μ[3,1]*Ω^-2 -c₂*(μ[1,0] + μ[1,1]) +
        c₃*μ[0,1]*Ω - c₄*μ[1,1] - c₁*μ[1,2]*Ω^-2 + c₁*μ[2,2]*Ω^-2
expr2 = simplify(value.(expr2))
@test isequal(expr1, expr2)

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
end
@test_throws ErrorException generate_raw_moment_eqs(rn, 2)

rn = @reaction_network begin
        (c₁*Y^3*X^1.5), X ⇒ Y
end
@test_throws ErrorException generate_raw_moment_eqs(rn, 2)

rn = @reaction_network begin
        (c₁*sin(exp(X)+Y)), X ⇒ Y
end
@test_throws ErrorException generate_raw_moment_eqs(rn, 2)

rn = @reaction_network begin
        (c₁*sin(exp(t))+Y^(3+1)), X ⇒ Y
end
eqs = generate_raw_moment_eqs(rn, 2)
@test eqs.q_order == 5