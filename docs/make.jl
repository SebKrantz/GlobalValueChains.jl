using GlobalValueChains
using Documenter

DocMeta.setdocmeta!(GlobalValueChains, :DocTestSetup, :(using GlobalValueChains); recursive = true)

makedocs(
    modules = [GlobalValueChains],
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
