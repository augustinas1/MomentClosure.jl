using MomentClosure, Symbolics, ModelingToolkit, Catalyst
using Test

@variables t, x(t)
@parameters κ, θ, σ

# Cox-Ingersoll-Ross model for interest rate forecasting 
# moment equations are naturally closed
cir_model = SDESystem([Differential(t)(x) ~ κ*(θ - x)], [σ * x^(1//2)], t, [x], [κ, θ, σ], name = :cir)
cir_moments = generate_raw_moment_eqs(cir_model, 2)
μ = cir_moments.μ
@test cir_moments.odes.eqs == [Differential(t)(μ[(1,)]) ~ expand(κ*(θ-μ[(1,)])), 
                               Differential(t)(μ[(2,)]) ~ expand(2*κ*(θ*μ[(1,)] - μ[(2,)]) + σ^2*μ[(1,)])]  

drift = [Differential(t)(x) ~ κ*(θ - x)]
diffusion = [σ * x^(1/2)]
ps = [κ, θ, σ]
cir_moments = generate_raw_moment_eqs(drift, diffusion, 2, ps)
μ = cir_moments.μ
@test cir_moments.odes.eqs == [Differential(t)(μ[(1,)]) ~ expand(κ*(θ-μ[(1,)])), 
                               Differential(t)(μ[(2,)]) ~ expand(2*κ*(θ*μ[(1,)] - μ[(2,)]) + σ^2*μ[(1,)])]  
   
cir_central_moments = generate_central_moment_eqs(cir_model, 2)
μ, M = cir_central_moments.μ, cir_central_moments.M
@test cir_central_moments.odes.eqs == [Differential(t)(μ[(1,)]) ~ expand(κ*(θ-μ[(1,)])),
                                       Differential(t)(M[(2,)]) ~ expand(σ^2*μ[(1,)] - 2*κ*M[(2,)])]

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
    μ = schloegl_moments.μ
    rhs(j) = j > 1 ? expand( j * ( k1*μ[(j-1,)] - k2*μ[(j,)] + (μ[(j+1,)] - μ[(j,)]) - (μ[(j+2,)] - 3*μ[(j+1,)] + 2*μ[(j,)]) ) 
                            + j*(j-1)/2 * ( k1*μ[(j-2,)] + k2*μ[(j-1,)] + (μ[(j,)] - μ[(j-1,)]) + (μ[(j+1,)] - 3*μ[(j,)] + 2*μ[(j-1,)]) ) ) : 
                    expand( k1*μ[(j-1,)] - k2*μ[(j,)] + (μ[(j+1,)] - μ[(j,)]) - (μ[(j+2,)] - 3*μ[(j+1,)] + 2*μ[(j,)]) ) 
    analytic_moment_eqs = [Differential(independent_variable(schloegl))(μ[(j,)]) ~ rhs(j) for j in 1:order]
    @test schloegl_moments.odes.eqs == analytic_moment_eqs
end

# check if things work for multiple species
rn = @reaction_network begin
    (c1), 2X + Y → 3X
    (c2), X → Y
    (c3, c4), 0 ↔ X
end c1 c2 c3 c4

c1, c2, c3, c4 = parameters(rn)

for order in 2:10
    rn_moments = generate_raw_moment_eqs(rn, order; langevin = true, combinatoric_ratelaw = false)
    μ = rn_moments.μ
    rhs(i,j) = (i >= 1 ? i*( c1*(μ[(i+1,j+1)]-μ[(i,j+1)]) - (c2 + c4)*μ[(i,j)] + c3*μ[(i-1,j)] ) : 0) +
            (j >= 1 ? j*( -c1*(μ[(i+2,j)]-μ[(i+1,j)]) + c2*μ[(i+1,j-1)] ) : 0) +
            1/2 * ( 
                ( i >= 2 ? i*(i-1) * ( c1*(μ[(i, j+1)] - μ[(i-1,j+1)]) + (c2 + c4)*μ[(i-1,j)] + c3*μ[(i-2,j)] ) : 0) +
                ( i >= 1 && j >= 1 ? 2*i*j * ( -c1*(μ[(i+1,j)] - μ[(i,j)]) - c2*μ[(i,j-1)] ) : 0 ) + 
                ( j >= 2 ? j*(j-1) * ( c1*(μ[(i+2,j-1)] - μ[(i+1,j-1)]) + c2*μ[(i+1,j-2)] ) : 0 )
            )
    analytic_moment_eqs = [Differential(independent_variable(rn))(μ[iter]) ~ expand(rhs(iter...)) for iter in vcat(rn_moments.iter_1, rn_moments.iter_m)]
    @test analytic_moment_eqs == rn_moments.odes.eqs
end
