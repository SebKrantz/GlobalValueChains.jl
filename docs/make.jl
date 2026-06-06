using ICIO
using Documenter

DocMeta.setdocmeta!(ICIO, :DocTestSetup, :(using ICIO); recursive = true)

makedocs(
    modules = [ICIO],
    authors = "Sebastian Krantz",
    sitename = "ICIO.jl",
    checkdocs = :exports,
    format = Documenter.HTML(
        canonical = "https://SebKrantz.github.io/ICIO.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "API reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/SebKrantz/ICIO.jl",
    devbranch = "main",
)
