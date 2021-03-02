using MomentClosure
using MomentClosure: value, define_M, define_μ, expand_mod
using Test
using Catalyst

@parameters t, k_on, k_off, k_p, γ_p, b
@variables p(t), g(t)

vars = [g, p]
ps = [k_on, k_off, k_p, γ_p, b]
S = [1 -1 0 0;
   0 0 b -1]
as = [k_on*(1-g),    # 0 -> g
    k_off*g*(p^2), # g -> 0
    k_p*g,         # 0 -> mP, m ~ Geometric(mean=b)
    γ_p*p]         # p -> 0
binary_vars = [1]
rn = ReactionSystemMod(t, vars, ps, as, S)

μ = define_μ(2,5)
M = define_M(2,5)

sys = generate_raw_moment_eqs(rn, 3)
expr1 = sys.odes.eqs[4].rhs
expr2 = k_on*μ[0,1] + b*k_p*μ[2,0] - k_off*μ[1,3] - k_on*μ[1,1]-γ_p*μ[1,1]
@test isequal(expr1, expr2)
@test length(sys.odes.eqs) == 9

sys_clean = bernoulli_moment_eqs(sys, binary_vars)
expr1 = sys_clean.odes.eqs[3].rhs
expr2 = k_on*μ[0,1] + b*k_p*μ[1,0] - k_off*μ[1,3] - k_on*μ[1,1]-γ_p*μ[1,1]
expr2 = simplify(value.(expr2))
@test isequal(expr1, expr2)
@test length(sys_clean.odes.eqs) == 6

@test_throws ErrorException moment_closure(sys, "conditional derivative matching", [])

closed_eqs = moment_closure(sys, "conditional gaussian", binary_vars)
expr1 = closed_eqs.closure[μ[1,4]]
expr2 = 4*μ[1,1]*μ[1,3]*μ[1,0]^-1 + 3*μ[1,2]^2*μ[1,0]^-1 - 12*μ[1,2]*μ[1,1]^2*μ[1,0]^-2 + 6*μ[1,1]^4*μ[1,0]^-3
@test isequal(expr1, expr2)
expr1 = closed_eqs.closure[μ[1,3]]
expr2 = 3*μ[1,2]*μ[1,1]*μ[1,0]^-1 - 2*μ[1,1]^3*μ[1,0]^-2
@test isequal(expr1, expr2)

closed_eqs = moment_closure(sys, "conditional derivative matching", binary_vars)
expr1 = closed_eqs.closure[μ[1,4]]
expr2 = μ[1,3]^4*μ[1,1]^4*μ[1,0]^-1*μ[1,2]^-6
@test isequal(expr1, expr2)
expr1 = closed_eqs.closure[μ[1,3]]
expr2 = μ[1,2]^3*μ[1,0]*μ[1,1]^-3
@test isequal(expr1, expr2)

# Testing formatting utilities
exprs = format_moment_eqs(closed_eqs)
# order of array elements can vary
ind = findfirst(!isnothing, match.(Ref(r"dμ₀₂/dt"), exprs))
expr1 = exprs[ind]
expr2 = "dμ₀₂/dt = γ_p*μ₀₁ + b*k_p*μ₁₀ + 2b*k_p*μ₁₁ + 2k_p*μ₁₀*(b^2) - (2γ_p*(μ₀₂))"
@test expr1 == expr2
exprs = format_closure(closed_eqs)
ind = findfirst(!isnothing, match.(Ref(r"μ₀₅ = "), exprs))
expr1 = exprs[ind]
expr2 = "μ₀₅ = ((μ₀₁)^-5)*((μ₀₂)^10)*((μ₀₃)^-10)*((μ₀₄)^5)"
@test expr1 == expr2

sys = generate_central_moment_eqs(rn, 3, 5)
expr1 = sys.odes.eqs[1].rhs
expr2 = -k_off*M[1,2] - k_on*μ[1,0] - k_off*M[0,2]*μ[1,0]-2*k_off*M[1,1]*μ[0,1] - k_off*μ[1,0]*μ[0,1]^2 + k_on
@test isequal(expand_mod(expr1), expr2)
@test length(sys.odes.eqs) == 9

sys_clean = bernoulli_moment_eqs(sys, binary_vars)
expr1 = expand_mod(sys_clean.odes.eqs[3].rhs)
expr2 = b*k_p*μ[1,0] - k_off*M[1,3] - k_on*M[1,1] - M[1,1]*γ_p - b*k_p*μ[1,0]^2 - k_off*M[0,3]*μ[1,0] -
      k_off*M[1,1]*μ[0,1]^2 - 2*k_off*M[1,2]*μ[0,1] - 2*k_off*M[0,2]*μ[0,1]*μ[1,0]
@test isequal(expr1, expr2)
@test length(sys_clean.odes.eqs) == 6

closed_eqs = moment_closure(sys, "conditional gaussian", binary_vars)
expr1 = expand_mod(closed_eqs.closure[M[0,5]])
expr2 = 45*μ[0,1]^5 + 10*M[0,2]*M[0,3] + 90*M[0,2]*μ[0,1]^3 + 5*M[0,4]*μ[0,1] + 30*μ[0,1]*M[0,2]^2 -
    45*μ[0,1]*(M[0,2]+μ[0,1]^2)^2
@test isequal(expr1, expr2)

closed_eqs = moment_closure(sys, "conditional derivative matching", binary_vars)
expr1 = expand_mod(closed_eqs.closure[M[1,3]])
expr2 = μ[1,0]*(M[1,1]+μ[0,1]*μ[1,0])^-3*(M[1,2]+M[0,2]*μ[1,0]+μ[1,0]*μ[0,1]^2+2*M[1,1]*μ[0,1])^3 -
    M[0,3]*μ[1,0] - 3*M[1,1]*μ[0,1]^2 - 3*M[1,2]*μ[0,1] - μ[1,0]*μ[0,1]^3- 3*M[0,2]*μ[0,1]*μ[1,0]
@test isequal(expr1, expr2)
