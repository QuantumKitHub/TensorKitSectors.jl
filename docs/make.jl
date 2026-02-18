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
        "Sector Interface" => [
            "Overview" => "interface/overview.md",
            "Required Methods" => "interface/required.md",
            "Optional Methods" => "interface/optional.md",
            "Traits and Styles" => "interface/traits.md",
            "Implementation Guidelines" => "interface/guidelines.md",
        ],
        "Sector Types" => [
            "Overview" => "sectors.md",
            "Abelian Groups" => [
                "Trivial" => "sectors/abelian/trivial.md",
                "ℤₙ (Cyclic)" => "sectors/abelian/zn.md",
                "U₁" => "sectors/abelian/u1.md",
            ],
            "Non-Abelian Groups" => [
                "SU₂" => "sectors/nonabelian/su2.md",
                "CU₁" => "sectors/nonabelian/cu1.md",
                "Dₙ (Dihedral)" => "sectors/nonabelian/dn.md",
                "A₄ (Alternating)" => "sectors/nonabelian/a4.md",
            ],
            "Anyonic Sectors" => [
                "PlanarTrivial" => "sectors/anyons/planartrivial.md",
                "Fibonacci" => "sectors/anyons/fibonacci.md",
                "Ising" => "sectors/anyons/ising.md",
            ],
            "Fermionic Sectors" => [
                "Fermion Parity" => "sectors/fermions/parity.md",
                "Fermion Number" => "sectors/fermions/number.md",
                "Fermion Spin" => "sectors/fermions/spin.md",
            ],
            "Composite Sectors" => [
                "Product" => "sectors/composite/product.md",
                "Time-Reversed" => "sectors/composite/timereversed.md",
            ],
        ],
        "Library" => "lib.md",
    ],
    checkdocs = :exports,
)

deploydocs(; repo = "github.com/QuantumKitHub/TensorKitSectors.jl.git", push_preview = true)
