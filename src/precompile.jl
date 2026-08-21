@setup_workload begin
    @compile_workload begin
        rn = Catalyst.@reaction_network begin
            @parameters k1 k2
            k1, 0 --> X
            k2, X --> 0
        end
        equations = generate_central_moment_eqs(rn, 2, 2)
        moment_closure(equations, "normal")
        generate_raw_moment_eqs(rn, 2)
    end
end
