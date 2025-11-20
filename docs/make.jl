using TensorKitSectors
using Documenter

mathengine = MathJax3(
    Dict(
        :loader => Dict("load" => ["[tex]/physics"]),
        :tex => Dict(
            "inlineMath" => [["\$", "\$"], ["\\(", "\\)"]],
            "tags" => "ams",
            "packages" => ["base", "ams", "autoload", "physics"]
        )
    )
)
makedocs(;
    sitename = "TensorKitSectors.jl",
    format = Documenter.HTML(;
        prettyurls = true,
        mathengine,
    ),
    pages = [
        "Home" => "index.md",
        "Library" => "lib.md",
    ],
    checkdocs = :exports,
)

deploydocs(; repo = "github.com/QuantumKitHub/TensorKitSectors.jl.git", push_preview = true)
