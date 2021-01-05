module CosmoCentral

using QuadGK
using Base: @kwdef
using Parameters
using NumericalIntegration
using Conda
Conda.add("numpy")
Conda.add("classy")
ENV["PYTHON"]=""
using Pkg
Pkg.build("PyCall")
using PyCall
numpy = pyimport("numpy")
classy = pyimport("classy")

include("CosmologicalParameters.jl")
include("Background.jl")
include("Density.jl")
include("Bias.jl")
include("MathUtils.jl")
include("WeightFunctions.jl")


export ComputeAdimensionalHubbleFactor

end # module
