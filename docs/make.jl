using Documenter, Corpuscles

makedocs(;
    modules = [Corpuscles],
    authors = "Tamas Gal and Johannes Schumann",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/Corpuscles_Logo.svg"],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ],
    repo="https://github.com/KM3NeT/Corpuscles.jl/blob/{commit}{path}#L{line}",
    sitename="Corpuscles.jl",
)

deploydocs(;
    repo="github.com/KM3NeT/Corpuscles.jl",
)
