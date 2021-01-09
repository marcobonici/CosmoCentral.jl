using Documenter, CosmoCentral
#include("/home/mbonici/Desktop/CosmoCentral.jl/src/CosmoCentral.jl")
using Plots, PlotThemes

# Set matplotlib gui backend
ENV["MPLBACKEND"] = "agg"

# Initialize backends
pyplot()
import PyPlot
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
font0 = Dict(
        "font.size" => 18,
        "axes.labelweight" => "bold",
        "axes.labelsize" => 16,
        "xtick.labelsize" => 8,
        "ytick.labelsize" => 32,
        "legend.fontsize" => 12,
)
merge!(rcParams, font0)

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
        "Cosmological Structure" => "CosmologicalStructure.md",
        "Math Utils" => "MathUtils.md"
    ]
)

deploydocs(
    repo = "github.com/marcobonici/CosmoCentral.jl.git",
    devbranch = "main"
)
