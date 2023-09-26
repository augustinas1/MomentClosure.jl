using MomentClosure, ModelingToolkit
using Documenter

makedocs(;
    modules=[MomentClosure, ModelingToolkit, Symbolics],
    authors="Augustinas Sukys",
    repo="https://github.com/augustinas1/MomentClosure.jl/blob/{commit}{path}#L{line}",
    sitename="MomentClosure.jl",
    doctest=false,
    warnonly = [:docs_block, :missing_docs, :cross_references],
    format=Documenter.HTML(;
        mathengine=Documenter.Writers.HTMLWriter.MathJax3(),
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://augustinas1.github.io/MomentClosure.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => Any[
            "theory/moment_expansion_CME.md",
            "theory/moment_expansion_SDE.md",
            "theory/moment_closure_approximations.md",
            "theory/linear_mapping_approximation.md"
        ],
        "Tutorials" => Any[
            "tutorials/using_momentclosure.md",
            "tutorials/using_momentclosure_SDE.md",
            "tutorials/common_issues.md",
            "tutorials/time-dependent_propensities.md",
            "tutorials/geometric_reactions+conditional_closures.md",
            "tutorials/P53_system_example.md",
            "tutorials/derivative_matching_example.md",
            "tutorials/SIR_example.md",
            "tutorials/LMA_example.md",
            "tutorials/parameter_estimation_SDE.md"
        ],
        "API" => Any[
            "api/momentclosure_api.md"
        ]
    ],

)

deploydocs(;
    repo="github.com/augustinas1/MomentClosure.jl",
    devbranch="main"
)
