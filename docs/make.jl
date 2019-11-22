using Documenter, Corpuscles

makedocs(;
    modules=[Corpuscles],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/KM3NeT/Corpuscles.jl/blob/{commit}{path}#L{line}",
    sitename="Corpuscles.jl",
    authors="Johannes Schumann, Tamas Gal",
    assets=String[],
)

deploydocs(;
    repo="github.com/KM3NeT/Corpuscles.jl",
)
