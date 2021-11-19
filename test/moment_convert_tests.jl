using MomentClosure
using MomentClosure: define_μ, define_M
#using SymbolicUtils: expand
using Symbolics: expand
using Test

M = define_M(2, 4)
κ = MomentClosure.cumulants_to_central_moments(2,4)
@test isequal(κ[1,2], M[1,2])
expr = M[3,1] - 3*M[1,1]*M[2,0]
@test isequal(κ[3,1], expr)

μ = define_μ(2, 3)
κ = MomentClosure.cumulants_to_raw_moments(2,3)
@test isequal(κ[1,0], μ[1,0])
expr = μ[2,1] + 2*μ[0,1]*μ[1,0]^2 - μ[0,1]*μ[2,0] - 2*μ[1,0]*μ[1,1]
# if this one fails most likely simplification functions have broke
@test isequal(expand(κ[2,1]),expr)

raw_to_central = MomentClosure.raw_to_central_moments(2,3)
expr = 2*μ[1,0]^3 + μ[3,0] - 3*μ[1,0]*μ[2,0]
@test isequal(raw_to_central[3,0], expr)
expr = μ[2,1] + 2*μ[0,1]*μ[1,0]^2 - 2*μ[1,0]*μ[1,1] - μ[0,1]*μ[2,0]
@test isequal(raw_to_central[2,1], expr)

central_to_raw = MomentClosure.central_to_raw_moments(2,3)
expr = M[2,0] + μ[1,0]^2
@test isequal(central_to_raw[2,0], expr)
expr = 2*M[1,1]*μ[0,1] + M[0,2]*μ[1,0] + M[1,2] + μ[1,0]*μ[0,1]^2
@test isequal(expand(central_to_raw[1,2]), expr)
