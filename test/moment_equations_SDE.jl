using MomentClosure, Symbolics, ModelingToolkit
using Test

@variables t, x(t)
@parameters κ, θ, σ

cir_model = SDESystem([Differential(t)(x) ~ κ*(θ - x)], [σ * x^(1/2)], t, [x], [κ, θ, σ], name = :cir)
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