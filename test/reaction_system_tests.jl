using MomentClosure
using MomentClosure: define_μ, define_M
using Test
using Catalyst

@parameters t, c₁, c₂, c₃, c₄, Ω
@variables X(t), Y(t)

S_mat = [ 1 -1  1 -1;
         -1  1  0  0]
# NOTE: only holds if combinatoric_ratelaw=false
a = [c₁*X*Y*(X-1)/Ω^2, c₂*X, c₃*Ω, c₄*X]

# checking if ReactionSystem and ReactionSystemMod are consistent

rn1 = @reaction_network begin
    (c₁/Ω^2), 2X + Y → 3X
    (c₂), X → Y
    (Ω*c₃, c₄), 0 ↔ X
end c₁ c₂ c₃ c₄ Ω

rn2 = ReactionSystemMod(t, [X, Y], [c₁, c₂, c₃, c₄, Ω], a, S_mat)

@test isequal(species(rn1), species(rn2))
@test params(rn1) == params(rn2)
@test speciesmap(rn1) == speciesmap(rn2)
@test paramsmap(rn1) == paramsmap(rn2)
@test numspecies(rn1) == numspecies(rn2) == 2
@test numreactions(rn1) == numreactions(rn2) == 4
@test numparams(rn1) == numparams(rn2) == 5
@test get_S_mat(rn1) == get_S_mat(rn2) == S_mat
@test isequal(propensities(rn1, combinatoric_ratelaw=false), propensities(rn2))
