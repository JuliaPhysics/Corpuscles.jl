using Documenter, Particles

makedocs(;
    modules=[Particles],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/KM3NeT/Particles.jl/blob/{commit}{path}#L{line}",
    sitename="Particles.jl",
    authors="Johannes Schumann, Tamas Gal",
    assets=String[],
)

deploydocs(;
    repo="github.com/KM3NeT/Particles.jl",
)
