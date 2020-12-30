using Documenter, CosmoCentral

push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "CosmoCentral.jl Documentation",
    authors="Marco Bonici",
    repo = "github.com/marcobonici/CosmoCentral.jl.git",
)
