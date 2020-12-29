using Documenter
using CosmoCentral

push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "CosmoCentral.jl Documentation",
    authors="Marco Bonici"
)
