using Documenter, Particles

makedocs(;
    modules=[Particles],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/8me/Particles.jl/blob/{commit}{path}#L{line}",
    sitename="Particles.jl",
    authors="Johannes Schumann, Tamas Gal",
    assets=String[],
)

deploydocs(;
    repo="github.com/8me/Particles.jl",
)
