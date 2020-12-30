using Documenter, CosmoCentral

push!(LOAD_PATH,"../src/")

makedocs(
    modules = [CosmoCentral],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "TaylorSeries.jl",
    authors  = "Luis Benet and David P. Sanders",
    pages = [
        "Home" => "index.md",
        "Background" => "anotherPage.md",
    ]
)

deploydocs(
    repo = "github.com/marcobonici/CosmoCentral.jl.git",
)
