using Documenter, CosmoCentral
using Plots, PlotThemes



push!(LOAD_PATH,"../src/")

makedocs(
    modules = [CosmoCentral],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename = "CosmoCentral.jl",
    authors  = "Marco Bonici",
    pages = [
        "Home" => "index.md",
        "Background Universe" => "BackgroundUniverse.md",
        "Source Density" => "SourceDensity.md",
        "Bias" => "Bias.md",
        "Boltzmann Solver" => "BoltzmannSolver.md",
        "Weight Functions" => "WeightFunction.md",
        "Angular Coefficients" => "AngularCoefficients.md",
        "Derivatives" => "Derivatives.md",
        "Cosmological Structure" => "CosmologicalStructure.md",
        "Math Utils" => "MathUtils.md"
    ]
)

deploydocs(
    repo = "github.com/marcobonici/CosmoCentral.jl.git",
    devbranch = "main"
)
