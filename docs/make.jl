using ICIO
using Documenter

DocMeta.setdocmeta!(ICIO, :DocTestSetup, :(using ICIO); recursive = true)

makedocs(
    modules = [ICIO],
    authors = "Sebastian Krantz",
    sitename = "GlobalValueChains.jl",
    checkdocs = :exports,
    format = Documenter.HTML(
        canonical = "https://SebKrantz.github.io/GlobalValueChains.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "API reference" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/SebKrantz/GlobalValueChains.jl",
    devbranch = "main",
)
