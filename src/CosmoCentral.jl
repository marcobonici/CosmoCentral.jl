__precompile__()
module CosmoCentral

using QuadGK: quadgk
using Base: @kwdef
using Parameters
using NumericalIntegration
using Interpolations
using PyCall
using Dierckx
using HDF5
using LoopVectorization
using Statistics
using JSON
using JSON3
using FFTW
using SpecialFunctions
using LatinHypercubeSampling
using Random
using DelimitedFiles
using LinearAlgebra
using Flux
numpy = pyimport("numpy")

const classy = PyNULL()


function __init__()
    copy!(classy, pyimport("classy"))
end

include("CosmologicalStructures.jl")
include("Background.jl")
include("Density.jl")
include("Bias.jl")
include("MathUtils.jl")
include("WeightFunctions.jl")
include("BoltzmannSolver.jl")
include("PowerSpectrum.jl")
include("AngularCoefficients.jl")
include("FileIO.jl")
include("Forecaster.jl")
include("Derivator.jl")
include("LHSampler.jl")
include("FFTLog.jl")
include("SourceFunctions.jl")
include("TransferFunctions.jl")
include("IntrinsicAlignment.jl")
include("Covariance.jl")
include("Fisher.jl")
include("MCMCUtils.jl")
include("Emulator.jl")
include("Probes.jl")

end # module
