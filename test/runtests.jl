using SafeTestsets

@safetestset "basic symbolics" begin include("basic_symbolic_tests.jl") end
@safetestset "moment conversion" begin include("moment_convert_tests.jl") end
@safetestset "reaction systems" begin include("reaction_system_tests.jl") end
@safetestset "moment eqs" begin include("moment_equations_tests.jl") end
@safetestset "moment closure" begin include("closure_tests.jl") end
@safetestset "moment closure (bernoulli)" begin include("bernoulli_closure_tests.jl") end
@safetestset "conditional closure" begin include("conditional_closure_tests.jl") end
@safetestset "latexify" begin include("latexify_tests.jl") end
@safetestset "sample moments" begin include("sample_moment_tests.jl") end
@safetestset "linear mapping approximation" begin include("lma_tests.jl") end
