using MomentClosure, Symbolics, ModelingToolkit, Catalyst
using ModelingToolkit: get_iv
using MomentClosure: define_μ, define_M, central_to_raw_moments
using Test

@variables t, x(t)
@parameters κ, θ, σ

# Cox-Ingersoll-Ross model for interest rate forecasting 
# moment equations are naturally closed
cir_model = SDESystem([Differential(t)(x) ~ κ*(θ - x)], [σ * x^(1//2)], t, [x], [κ, θ, σ], name = :cir)

cir_moments = generate_raw_moment_eqs(cir_model, 2)
μ = cir_moments.μ
@test isequal(cir_moments.odes.eqs,  [Differential(t)(μ[(1,)]) ~ expand(κ*(θ-μ[(1,)])), 
                                      Differential(t)(μ[(2,)]) ~ expand(2*κ*(θ*μ[(1,)] - μ[(2,)]) + σ^2*μ[(1,)])])

cir_central_moments = generate_central_moment_eqs(cir_model, 2)
@test isequal(cir_central_moments.odes, generate_central_moment_eqs(cir_model, 2, 2).odes)
μ, M = cir_central_moments.μ, cir_central_moments.M
@test isequal(cir_central_moments.odes.eqs, [Differential(t)(μ[(1,)]) ~ expand(κ*(θ-μ[(1,)])),
                                             Differential(t)(M[(2,)]) ~ expand(σ^2*μ[(1,)] - 2*κ*M[(2,)])])

## Chemical Reaction Networks via Chemical Langevin Equation
# unimolecular system
schloegl = @reaction_network begin
    1.0, 2X → 3X
    1.0, 3X → 2X
    k1, ∅ → X
    k2, X → ∅
end k1 k2

k1, k2 = parameters(schloegl)
for order in 2:6
    schloegl_moments = generate_raw_moment_eqs(schloegl, order; langevin = true, combinatoric_ratelaw = false)
    local μ = schloegl_moments.μ
    rhs(j) = j > 1 ? expand( j * ( k1*μ[(j-1,)] - k2*μ[(j,)] + (μ[(j+1,)] - μ[(j,)]) - (μ[(j+2,)] - 3*μ[(j+1,)] + 2*μ[(j,)]) ) 
                            + j*(j-1)/2 * ( k1*μ[(j-2,)] + k2*μ[(j-1,)] + (μ[(j,)] - μ[(j-1,)]) + (μ[(j+1,)] - 3*μ[(j,)] + 2*μ[(j-1,)]) ) ) : 
                    expand( k1*μ[(j-1,)] - k2*μ[(j,)] + (μ[(j+1,)] - μ[(j,)]) - (μ[(j+2,)] - 3*μ[(j+1,)] + 2*μ[(j,)]) ) 
    analytic_moment_eqs = [Differential(get_iv(schloegl))(μ[(j,)]) ~ rhs(j) for j in 1:order]
    @test isequal(schloegl_moments.odes.eqs, analytic_moment_eqs)
end

μ = define_μ(1, 4); M = define_M(1, 4)
schloegl_moments = generate_central_moment_eqs(schloegl, 2, 4, langevin=true, combinatoric_ratelaw = false)
expr = k1 + 4*μ[(1,)]^2 + 4*M[(2,)] - μ[(1,)]^3 - M[(3,)] - 3*μ[(1,)] - k2*μ[(1,)] - 3*M[(2,)]*μ[(1,)]
@test isequal(schloegl_moments.odes.eqs[1].rhs, expr)
expr = k1 + 9*M[(3,)] + k2*μ[(1,)] + μ[(1,)]^3 + 19*M[(2,)]*μ[(1,)] + μ[(1,)] - 2*μ[(1,)]^2 - 
       2*M[(4,)] - 8*M[(2,)] - 2*k2*M[(2,)] - 6*μ[(1,)]^2*M[(2,)] - 6*M[(3,)]*μ[(1,)]
@test isequal(schloegl_moments.odes.eqs[2].rhs, expr)

# check if things work for multiple species
rn = @reaction_network begin
    (c1), 2X + Y → 3X
    (c2), X → Y
    (c3, c4), 0 ↔ X
end c1 c2 c3 c4

c1, c2, c3, c4 = parameters(rn)

for order in 2:10
    rn_moments = generate_raw_moment_eqs(rn, order; langevin = true, combinatoric_ratelaw = false)
    local μ = rn_moments.μ
    rhs(i,j) = (i >= 1 ? i*( c1*(μ[(i+1,j+1)]-μ[(i,j+1)]) - (c2 + c4)*μ[(i,j)] + c3*μ[(i-1,j)] ) : 0) +
            (j >= 1 ? j*( -c1*(μ[(i+2,j)]-μ[(i+1,j)]) + c2*μ[(i+1,j-1)] ) : 0) +
            1/2 * ( 
                ( i >= 2 ? i*(i-1) * ( c1*(μ[(i, j+1)] - μ[(i-1,j+1)]) + (c2 + c4)*μ[(i-1,j)] + c3*μ[(i-2,j)] ) : 0) +
                ( i >= 1 && j >= 1 ? 2*i*j * ( -c1*(μ[(i+1,j)] - μ[(i,j)]) - c2*μ[(i,j-1)] ) : 0 ) + 
                ( j >= 2 ? j*(j-1) * ( c1*(μ[(i+2,j-1)] - μ[(i+1,j-1)]) + c2*μ[(i+1,j-2)] ) : 0 )
            )
    analytic_moment_eqs = [Differential(get_iv(rn))(μ[iter]) ~ expand(rhs(iter...)) for iter in vcat(rn_moments.iter_1, rn_moments.iter_m)]
    isequal(analytic_moment_eqs[1], rn_moments.odes.eqs[1])
    @test isequal(analytic_moment_eqs, rn_moments.odes.eqs)
end

raw_eqs = generate_raw_moment_eqs(rn, 2, langevin = true, combinatoric_ratelaw=false)
central_eqs = generate_central_moment_eqs(rn, 2, langevin=true, combinatoric_ratelaw=false)
μ = define_μ(2, 4); M = define_μ(2, 4)
central_to_raw = central_to_raw_moments(2, 4)
subdict = Dict( μ[iter] => central_to_raw[iter] for iter in filter(x -> sum(x) > 1, raw_eqs.iter_all))

expr = expand(substitute(raw_eqs.odes.eqs[1].rhs, subdict))
@test isequal(expr, central_eqs.odes.eqs[1].rhs)
expr = expand(substitute( raw_eqs.odes.eqs[5].rhs - 2*μ[0,1]*raw_eqs.odes.eqs[2].rhs, subdict))
@test isequal(expr, central_eqs.odes.eqs[5].rhs)
expr = expand(substitute( raw_eqs.odes.eqs[4].rhs - μ[1,0]*raw_eqs.odes.eqs[2].rhs - μ[0,1]*raw_eqs.odes.eqs[1].rhs, subdict ))
@test isequal(expr, central_eqs.odes.eqs[4].rhs)

#=
using OrdinaryDiffEq
using ModelingToolkit: modelingtoolkitize
f(du,u,p,t) = du .= 1.01u
function g(du,u,p,t)
  du[1,1] = 0.3u[1]
  du[1,2] = 0.6u[1]
  du[1,3] = 0.9u[1]
  du[1,4] = 0.12u[1]
  du[2,1] = 1.2u[2]
  du[2,2] = 0.2u[2]
  du[2,3] = 0.3u[2]
  du[2,4] = 1.8u[2]
end
prob = SDEProblem(f,g,ones(2),(0.0,1.0),noise_rate_prototype=zeros(2,4)) 
sys = modelingtoolkitize(prob)
eqs = generate_raw_moment_eqs(sys, 2)
=#