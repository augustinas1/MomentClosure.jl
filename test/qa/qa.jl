using SciMLTesting, MomentClosure, Test
using JET

run_qa(
    MomentClosure;
    explicit_imports = true,
    ei_kwargs = (;
        # Genuine non-public internals of non-SciML / not-yet-`public`-declared
        # deps that MomentClosure explicitly imports. Re-checked empirically on
        # Julia 1.12 against the registered make-public releases (SymbolicUtils
        # 4.37, Symbolics 7.29, TupleTools 1.6): each of these is still flagged
        # as non-public, and none has a public re-exporting owner to migrate to.
        all_explicit_imports_are_public = (;
            ignore = (
                :FnType, :ispow,                         # SymbolicUtils internals
                :map_subscripts, :var_from_nested_derivative,  # Symbolics internals
                :sort,                                   # TupleTools (non-SciML dep)
            ),
        ),
    ),
    # Resolving the ~45 implicit imports would require enumerating hundreds of
    # Catalyst/ModelingToolkit/Symbolics re-exports by hand, which is fragile
    # across the MTK 11 / Symbolics 7 line. Tracked instead of mass-refactored.
    # https://github.com/SciML/MomentClosure.jl/issues/101
    ei_broken = (:no_implicit_imports,),
)
