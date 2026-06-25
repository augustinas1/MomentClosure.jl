using SciMLTesting, MomentClosure, Test
using JET

run_qa(
    MomentClosure;
    explicit_imports = true,
    ei_kwargs = (;
        # Names imported/accessed from a re-exporting module rather than their
        # defining owner. These are other packages' (re)exports whose owner moved
        # as ModelingToolkit/Symbolics/SymbolicUtils refactored; MomentClosure
        # uses the public-facing re-export, so the owner mismatch is expected.
        all_explicit_imports_via_owners = (;
            ignore = (
                :get_eqs, :get_iv, :get_noise_eqs, :get_ps,  # owner ModelingToolkitBase, re-exported by ModelingToolkit
                :getname,                                     # owner SymbolicIndexingInterface, re-exported by ModelingToolkit
                :scalarize,                                   # owner SymbolicUtils, re-exported by Symbolics
            ),
        ),
        all_qualified_accesses_via_owners = (;
            ignore = (
                :get_eqs, :get_iv, :get_ps,  # owner ModelingToolkitBase, accessed via ModelingToolkit
            ),
        ),
        # Other packages' non-public (not yet `public`-declared) names that
        # MomentClosure genuinely depends on; these become public as the base
        # libraries (ModelingToolkit/Symbolics/SymbolicUtils/SciMLBase/TupleTools)
        # tag releases that mark them.
        all_explicit_imports_are_public = (;
            ignore = (
                :BasicSymbolic, :FnType, :isadd, :isdiv, :ismul, :ispow, :isterm,  # SymbolicUtils
                :map_subscripts, :scalarize, :value, :var_from_nested_derivative,  # Symbolics
                :get_eqs, :get_iv, :get_noise_eqs, :get_ps, :getname,              # ModelingToolkit
                :NullParameters,                                                   # SciMLBase
                :sort,                                                             # TupleTools
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                :VariableSource, :variable,  # Symbolics
                :get_eqs, :get_iv, :get_ps,  # ModelingToolkit
                :filter,                     # Base.Iterators
            ),
        ),
    ),
    # Resolving the ~45 implicit imports would require enumerating hundreds of
    # Catalyst/ModelingToolkit/Symbolics re-exports by hand, which is fragile
    # across the MTK 11 / Symbolics 7 line. Tracked instead of mass-refactored.
    # https://github.com/SciML/MomentClosure.jl/issues/101
    ei_broken = (:no_implicit_imports,),
)
