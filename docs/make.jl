using Documenter, CosmoCentral
#include("/home/mbonici/Desktop/CosmoCentral.jl/src/CosmoCentral.jl")

push!(LOAD_PATH,"../src/")

makedocs(
    modules = [CosmoCentral],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "CosmoCentral.jl",
    authors  = "Marco Bonici",
    pages = [
        "Home" => "index.md",
        "Background" => "anotherPage.md",
    ]
)

deploydocs(
    repo = "github.com/marcobonici/CosmoCentral.jl.git",
    devbranch = "main"
)
