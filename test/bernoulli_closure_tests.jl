using MomentClosure
using MomentClosure: define_M, define_μ
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

μ = define_μ(2,4)
M = define_M(2,4)
sys = generate_central_moment_eqs(rn, 2)

closed_eqs = moment_closure(sys, "zero", binary_vars)
@test length(closed_eqs.odes.states) == 4 && isequal(closed_eqs.closure[M[0,4]], 0)

closed_eqs = moment_closure(sys, "normal", binary_vars)
@test length(closed_eqs.odes.states) == 4 && isequal(closed_eqs.closure[M[0,4]], 3*M[0,2]^2)

closed_eqs= moment_closure(sys, "log-normal", binary_vars)
expr1 = closed_eqs.closure[M[0,3]]
expr2 = μ[0,1]^3*(1+M[0,2]*μ[0,1]^-2)^3 - μ[0,1]^3 - 3*M[0,2]*μ[0,1]
@test length(closed_eqs.odes.states) == 4 && isequal(expand(expr1), expr2)

closed_eqs = moment_closure(sys, "poisson", binary_vars)
@test length(closed_eqs.odes.states) == 4 && isequal(closed_eqs.closure[M[0,3]], μ[0,1])

closed_eqs = moment_closure(sys, "gamma", binary_vars)
expr1 = closed_eqs.closure[M[0,3]]
expr2 = 2*M[0,2]^2*μ[0,1]^-1
@test length(closed_eqs.odes.states) == 4 && isequal(expr1, expr2)

closed_eqs = moment_closure(sys, "derivative matching", binary_vars)
expr1 = closed_eqs.closure[sys.M[0,3]]
expr2 = μ[0,1]^-3 * (M[0,2]+μ[0,1]^2)^3 - μ[0,1]^3 - 3*M[0,2]*μ[0,1]
@test length(closed_eqs.odes.states) == 4 && isequal(expr1, expr2)

sys = generate_raw_moment_eqs(rn, 2)

closed_eqs = moment_closure(sys, "zero", binary_vars)
@test length(closed_eqs.odes.states) == 4

closed_eqs = moment_closure(sys, "normal", binary_vars)
@test length(closed_eqs.odes.states) == 4

closed_eqs = moment_closure(sys, "log-normal", binary_vars)
@test length(closed_eqs.odes.states) == 4

closed_eqs = moment_closure(sys, "poisson", binary_vars)
@test length(closed_eqs.odes.states) == 4

closed_eqs = moment_closure(sys, "gamma", binary_vars)
@test isequal(expr1, expr2)

closed_eqs = moment_closure(sys, "derivative matching", binary_vars)
@test isequal(expr1, expr2)
