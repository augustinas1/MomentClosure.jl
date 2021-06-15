using MomentClosure, ModelingToolkit
using Documenter

makedocs(;
    modules=[MomentClosure, ModelingToolkit, Symbolics],
    authors="Augustinas Sukys",
    repo="https://github.com/augustinas1/MomentClosure.jl/blob/{commit}{path}#L{line}",
    sitename="MomentClosure.jl",
    doctest=false,
    format=Documenter.HTML(;
        mathengine=Documenter.Writers.HTMLWriter.MathJax2(),
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://augustinas1.github.io/MomentClosure.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => Any[
            "theory/moment_expansion.md",
            "theory/moment_closure_approximations.md"
        ],
        "Tutorials" => Any[
            "tutorials/using_momentclosure.md",
            "tutorials/common_issues.md",
            "tutorials/time-dependent_propensities.md",
            "tutorials/geometric_reactions+conditional_closures.md",
            "tutorials/P53_system_example.md",
            "tutorials/derivative_matching_example.md"
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
