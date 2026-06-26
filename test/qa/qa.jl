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
        # Other packages' still-non-public names that MomentClosure imports.
        # These remain after the base libs' public-API releases (SymbolicUtils
        # 4.36, Symbolics 7.28.1, ModelingToolkit 11.29, SciMLBase 3.24): they
        # are genuine internals not yet `public`-declared in the importing module
        # (plus `sort` from TupleTools, a non-SciML dep). Re-checked on Julia 1.12.
        all_explicit_imports_are_public = (;
            ignore = (
                :FnType, :isdiv, :ispow, :isterm,        # SymbolicUtils
                :map_subscripts, :var_from_nested_derivative,  # Symbolics
                :getname,                                # ModelingToolkit (owner SymbolicIndexingInterface)
                :sort,                                   # TupleTools
            ),
        ),
        all_qualified_accesses_are_public = (; ignore = ()),
    ),
    # Resolving the ~45 implicit imports would require enumerating hundreds of
    # Catalyst/ModelingToolkit/Symbolics re-exports by hand, which is fragile
    # across the MTK 11 / Symbolics 7 line. Tracked instead of mass-refactored.
    # https://github.com/SciML/MomentClosure.jl/issues/101
    ei_broken = (:no_implicit_imports,),
)
