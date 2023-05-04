using Test
using MomentClosure
using Catalyst

@parameters t, c₁, c₂, c₃, c₄, Ω
@species X(t), Y(t)

S_mat = [ 1 -1  1 -1;
         -1  1  0  0]
# NOTE: only holds if combinatoric_ratelaw=false
a = [c₁*X*Y*(X-1)/Ω^2, c₂*X, c₃*Ω, c₄*X]

rn = @reaction_network begin
    @parameters c₁ c₂ c₃ c₄ Ω
    (c₁/Ω^2), 2X + Y → 3X
    (c₂), X → Y
    (Ω*c₃, c₄), 0 ↔ X
end
smap = speciesmap(rn)

@test isequal(get_stoichiometry(rn, smap), S_mat)
@test isequal(propensities(rn, combinatoric_ratelaw=false), a)
