using Documenter
using CosmoCentral

push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "CosmoCentral.jl Documentation",
    pages = [
       "Index" => "index.md",
       "anotherPage" => "anotherPage.md"
    ],
    format = Documenter.HTML(),
    modules = [CosmoCentral]
)

deploydocs(
    repo = "github.com/marcobonici/CosmoCentral.jl.git",
    devbranch = "main"
)
