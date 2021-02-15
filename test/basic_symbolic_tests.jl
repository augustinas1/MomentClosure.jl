using MomentClosure
using MomentClosure: define_μ, define_M
using Test, SafeTestsets

@testset "creating symbolic moments" begin

    μ = define_μ(3, 3)
    @test μ[(0,0,0)] == 1
    @test string(μ[(2,1,0)]) == "μ₂₁₀(t)"
    @test string(μ[(0,0,3)]) == "μ₀₀₃(t)"
    M = define_M(3, 3)
    @test M[(1,0,0)] == 0
    @test string(M[1,0,2]) == "M₁₀₂(t)"

end
