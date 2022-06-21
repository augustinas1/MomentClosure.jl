using MomentClosure
using SciMLBase, JLD2
using Test

JLD2.@load joinpath(@__DIR__, "test_data.jld2") sol

μ_t = Dict([(0, 1, 0) => [30.0, 33.6, 37.6]
            (0, 0, 2) => [3600.0, 3461.4, 3207.0]
            (0, 0, 1) => [60.0, 58.8, 56.6]
            (1, 0, 0) => [70.0, 68.4, 71.0]
            (0, 2, 0) => [900.0, 1134.2, 1420.8]
            (1, 0, 1) => [4200.0, 4021.8, 4018.9]
            (0, 1, 1) => [1800.0, 1973.8, 2124.2]
            (2, 0, 0) => [4900.0, 4687.2, 5058.2]
            (1, 1, 0) => [2100.0, 2299.1, 2668.3]])

μ_1 = get_raw_moments(sol, 2)
μ_2 = get_raw_moments(sol, 2, naive=false, b=3)
x_1 = [μ_t[key] ≈ μ_1[key] for key in keys(μ_t)]
x_2 = [μ_t[key] ≈ μ_2[key] for key in keys(μ_t)]
@test all(x_1) && all(x_2)

M_t = Dict([(0, 0, 2) => [0.0, 3.96, 3.44]
            (0, 1, 1) => [0.0, -1.88, -3.96]
            (0, 2, 0) => [0.0, 5.24, 7.04]
            (1, 0, 1) => [0.0, -0.12, 0.3]
            (2, 0, 0) => [0.0, 8.64, 17.2]
            (1, 1, 0) => [0.0, 0.86, -1.3]])

M_1 = get_central_moments(sol, 2)
M_2 = get_central_moments(sol, 2, naive=false, b=3)
x_1 = [M_t[key] ≈ M_1[key] for key in keys(M_t)]
x_2 = [M_t[key] ≈ M_2[key] for key in keys(M_t)]
@test all(x_1) && all(x_2)

κ_t = Dict([(0, 1, 0) => [30.0, 33.6, 37.6]
            (0, 0, 2) => [0.0, 3.96, 3.44]
            (0, 0, 1) => [60.0, 58.8, 56.6]
            (1, 0, 0) => [70.0, 68.4, 71.0]
            (0, 2, 0) => [0.0, 5.24, 7.04]
            (1, 0, 1) => [0.0, -0.12, 0.3]
            (0, 1, 1) => [0.0, -1.88, -3.96]
            (2, 0, 0) => [0.0, 8.64, 17.2]
            (1, 1, 0) => [0.0, 0.86, -1.3]])

κ_1 = get_cumulants(sol, 2)
κ_2 = get_cumulants(sol, 2, naive=false, b=3)
x_1 = [κ_t[key] ≈ κ_1[key] for key in keys(κ_t)]
x_2 = [κ_t[key] ≈ κ_2[key] for key in keys(κ_t)]
@test all(x_1) && all(x_2)