using MomentClosure
using MomentClosure: define_M, define_μ
using Test
using Catalyst
using Distributions
using Symbolics: substitute, value, get_variables

# Recent Symbolics keeps expressions such as `x^-1` and `(1 / x)` in distinct (but
# mathematically equal) forms, so `isequal`/`simplify` no longer reliably proves two
# closure expressions equal. Check equality by evaluating the difference at a generic
# point instead; identically-zero rational functions vanish at any non-singular point.
function symeq(expr1, expr2; atol = 1.0e-10)
    diff = expr1 - expr2
    vars = get_variables(diff)
    point = Dict(v => 1.0 + i / 17 for (i, v) in enumerate(vars))
    return isapprox(Float64(value(substitute(diff, point))), 0; atol = atol)
end

@parameters p
@register_symbolic Distributions.Geometric(p)
m = rand(Distributions.Geometric(p))

rn = @reaction_network begin
    k_on * (1 - g), 0 --> g
    k_off * P^2, g --> 0
    k_p, g --> g + $m * P
    γ_p, P --> 0
end

binary_vars = [1]

μ = define_μ(2, 4)
M = define_M(2, 4)
sys = generate_central_moment_eqs(rn, 2, combinatoric_ratelaws = false)

closed_eqs = moment_closure(sys, "zero", binary_vars)
@test length(unknowns(closed_eqs)) == 4 && isequal(get_closure(closed_eqs)[M[0, 4]], 0)

closed_eqs = moment_closure(sys, "normal", binary_vars)
@test length(unknowns(closed_eqs)) == 4 && isequal(get_closure(closed_eqs)[M[0, 4]], 3 * M[0, 2]^2)

closed_eqs = moment_closure(sys, "log-normal", binary_vars)
expr1 = expand(get_closure(closed_eqs)[M[0, 3]])
expr2 = M[0, 2]^3 * μ[0, 1]^-3 + 3 * M[0, 2]^2 * μ[0, 1]^-1
@test length(unknowns(closed_eqs)) == 4 && iszero(simplify(expr1 - expr2))

closed_eqs = moment_closure(sys, "poisson", binary_vars)
@test length(unknowns(closed_eqs)) == 4 && isequal(get_closure(closed_eqs)[M[0, 3]], μ[0, 1])

closed_eqs = moment_closure(sys, "gamma", binary_vars)
expr1 = get_closure(closed_eqs)[M[0, 3]]
expr2 = 2 * M[0, 2]^2 * μ[0, 1]^-1
@test length(unknowns(closed_eqs)) == 4 && symeq(expr1, expr2)

closed_eqs = moment_closure(sys, "derivative matching", binary_vars)
expr1 = expand(get_closure(closed_eqs)[sys.M[0, 3]])
expr2 = M[0, 2]^3 * μ[0, 1]^-3 + 3 * M[0, 2]^2 * μ[0, 1]^-1
@test length(unknowns(closed_eqs)) == 4 && iszero(simplify(expr1 - expr2))

sys = generate_raw_moment_eqs(rn, 2)

closed_eqs = moment_closure(sys, "zero", binary_vars)
@test length(unknowns(closed_eqs)) == 4

closed_eqs = moment_closure(sys, "normal", binary_vars)
@test length(unknowns(closed_eqs)) == 4

closed_eqs = moment_closure(sys, "log-normal", binary_vars)
@test length(unknowns(closed_eqs)) == 4

closed_eqs = moment_closure(sys, "poisson", binary_vars)
@test length(unknowns(closed_eqs)) == 4

closed_eqs = moment_closure(sys, "gamma", binary_vars)
@test length(unknowns(closed_eqs)) == 4

closed_eqs = moment_closure(sys, "derivative matching", binary_vars)
@test length(unknowns(closed_eqs)) == 4
