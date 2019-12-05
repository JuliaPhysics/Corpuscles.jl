using Documenter, Corpuscles

makedocs(;
    modules = [Corpuscles],
    authors = "Tamas Gal and Johannes Schumann",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/logo.svg"],
    ),
    pages=[
        "Introduction" => "readme.md",
        "Acknowledgements" => "acknowledgements.md",
        "API" => "api.md",
    ],
    repo="https://github.com/KM3NeT/Corpuscles.jl/blob/{commit}{path}#L{line}",
    sitename="Corpuscles.jl",
)

deploydocs(;
    repo="github.com/KM3NeT/Corpuscles.jl",
)
