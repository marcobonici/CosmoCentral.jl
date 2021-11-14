using Documenter
using CosmoCentral
using Plots
using PlotThemes

ENV["GKSwstype"] = "100"

push!(LOAD_PATH,"../src/")

makedocs(
    modules = [CosmoCentral],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true",
    sidebar_sitename=false),
    sitename = "CosmoCentral.jl",
    authors  = "Marco Bonici",
    pages = [
        "Home" => "index.md",
        "Cosmology" => "Cosmology.md",
        "Background Universe" => "BackgroundUniverse.md",
        "Source Density" => "SourceDensity.md",
        "Bias" => "Bias.md",
        "Boltzmann Solver" => "BoltzmannSolver.md",
        "Weight Functions" => "WeightFunction.md",
        "Angular Coefficients" => "AngularCoefficients.md",
        "Derivatives" => "Derivatives.md",
        "Covariance Matrix" => "Covariance.md",
        "Fisher Matrix" => "Fisher.md",
        "Math Utils" => "MathUtils.md"
    ]
)

deploydocs(
    repo = "github.com/marcobonici/CosmoCentral.jl.git",
    devbranch = "develop"
)
