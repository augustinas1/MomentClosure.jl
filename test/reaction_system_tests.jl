using MomentClosure
using MomentClosure: define_μ, define_M, isvar
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
ps1 = reactionparams(rn1); ps2 = reactionparams(rn2)
@test length(ps1) == length(ps2) == 5
@test all(isvar.(ps1, Ref(ps2)))
@test all(isvar.(ps2, Ref(collect(keys(paramsmap(rn2))))))
@test isequal(speciesmap(rn1), speciesmap(rn2))
@test numspecies(rn1) == numspecies(rn2) == 2
@test numreactions(rn1) == numreactions(rn2) == 4
@test netstoichmat(rn1) == netstoichmat(rn2) == S_mat
@test isequal(propensities(rn1, combinatoric_ratelaw=false), propensities(rn2))
