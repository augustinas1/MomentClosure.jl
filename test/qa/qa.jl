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
        # Genuine non-public internals of non-SciML / not-yet-`public`-declared
        # deps that MomentClosure explicitly imports. Re-checked empirically on
        # Julia 1.12 against the registered make-public releases (SymbolicUtils
        # 4.37, Symbolics 7.29, ModelingToolkit 11.29): each of these is still
        # flagged as non-public.
        all_explicit_imports_are_public = (;
            ignore = (
                :FnType, :ispow,                         # SymbolicUtils internals
                :map_subscripts, :var_from_nested_derivative,  # Symbolics internals
                :sort,                                   # TupleTools (non-SciML dep)
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
