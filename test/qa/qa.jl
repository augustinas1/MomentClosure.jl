using MomentClosure
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(MomentClosure)
end

@testset "JET" begin
    JET.test_package(MomentClosure; target_defined_modules = true)
end
