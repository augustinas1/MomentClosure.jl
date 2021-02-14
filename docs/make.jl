using MomentClosure
using Documenter

makedocs(;
    modules=[MomentClosure],
    authors="Augustinas Sukys",
    repo="https://github.com/augustinas1/MomentClosure.jl/blob/{commit}{path}#L{line}",
    sitename="MomentClosure.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://augustinas1.github.io/MomentClosure.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/augustinas1/MomentClosure.jl",
)
