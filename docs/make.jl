using MomentClosure
using Documenter

makedocs(;
    modules=[MomentClosure],
    authors="Augustinas Sukys",
    repo="https://github.com/augustinas1/MomentClosure.jl/blob/{commit}{path}#L{line}",
    sitename="MomentClosure.jl",
    format=Documenter.HTML(;
        mathengine=Documenter.Writers.HTMLWriter.MathJax2(),
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://augustinas1.github.io/MomentClosure.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => Any[
            "theory/moment_equations.md",
            "theory/moment_closure_approximations.md"
        ],
        "Tutorials" => Any[
            "tutorials/model_initialisation.md",
            "tutorials/SSA.md",
            "tutorials/moment_closure.md",
            "tutorials/Brusselator.md",
            "tutorials/P53_system.md",
            "tutorials/derivative_matching_example.md",
            "tutorials/gene_network.md"
        ]
    ],
)

deploydocs(;
    repo="github.com/augustinas1/MomentClosure.jl",
    devbranch="main"
)
