using MomentClosure
using MomentClosure: define_μ
using Symbolics: value
using Catalyst
using Test

@parameters t, σ_b, σ_u, ρ_b, ρ_u, σ_b_LMA
@variables g(t)

# simple feedback loop
rn_nonlinear = @reaction_network begin
      σ_b, g + p → 0
      σ_u*(1-g), 0 ⇒ g + p
      ρ_u, g → g + p
      ρ_b*(1-g), 0 ⇒ p
      1, p → 0
end

rn_linear = @reaction_network begin
      σ_b_LMA, g → 0
      σ_u*(1-g), 0 ⇒ g
      ρ_u, g → g+p
      (ρ_b*(1-g)), 0 ⇒ p
      1, p → 0
end

binary_vars = [speciesmap(rn_nonlinear)[g]]
LMA_eqs, _ = linear_mapping_approximation(rn_nonlinear, rn_linear, binary_vars, combinatoric_ratelaws=false)

μ = define_μ(2,3)
expr1 = LMA_eqs.odes.eqs[1].rhs
expr2 = σ_u - σ_b*μ[1,1] - σ_u*μ[1,0]
@test isequal(expr1, value.(expr2))
expr1 = LMA_eqs.odes.eqs[2].rhs
expr2 = ρ_b + ρ_u*μ[1,0] - ρ_b*μ[1,0] - μ[0,1]
@test isequal(expr1, value.(expr2))
expr1 = LMA_eqs.odes.eqs[3].rhs
expr2 = ρ_u*μ[1,0] + σ_u*μ[0,1] - μ[1,1] - σ_u*μ[1,1] - σ_b*μ[1,1]^2*μ[1,0]^-1
@test isequal(simplify(expr1), simplify(expr2))

# cooperativity cp=2
rn_nonlinear = @reaction_network begin
      σ_b, g + 2*p → 0
      σ_u*(1-g), 0 ⇒ g + 2*p
      ρ_u, g → g + p
      ρ_b*(1-g), 0 ⇒ p
      1, p → 0
end

binary_vars = [speciesmap(rn_nonlinear)[g]]
_, effective_params = linear_mapping_approximation(rn_nonlinear, rn_linear, binary_vars, combinatoric_ratelaws=false)
expr1 = effective_params[value(σ_b_LMA)]
expr2 = (σ_b*μ[1,2] - σ_b*μ[1,1]) * μ[1,0]^-1
@test isequal(expr1, value(expr2))