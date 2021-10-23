"""
    ComputeBias(z::Float64, ::PiecewiseBias,
    ConvolvedDensity::AbstractConvolvedDensity)

This function evaluate the piecewise bias, given by ``\\sqrt{1+z}``, for a given
redshift.

"""
function ComputeBias(z::T, Piecewise::PiecewiseBias,
    ConvolvedDensity::AbstractConvolvedDensity) where T
    idx = BinSearch(z, ConvolvedDensity.ZBinArray)
    bias = sqrt(1+(ConvolvedDensity.ZBinArray[idx]+
    ConvolvedDensity.ZBinArray[idx+1])/2) * Piecewise.BiasMultiplier[idx]
    return bias
end

function ComputeBias(z::T, Bias::EuclidBias,
    ConvolvedDensity::AbstractConvolvedDensity) where T
    bias = Bias.A+Bias.B/(1+exp((Bias.D-z)*Bias.C))
    return bias
end

"""
    ComputeBiasGrid!(CosmologicalGrid::CosmologicalGrid,
    gcWeightFunction::GCWeightFunction, Bias::AbstractBias,
    ConvolvedDensity::AbstractConvolvedDensity)

This function evaluate the bias, over the cosmological redshift grid.
"""
function ComputeBiasGrid!(cosmologicalGrid::CosmologicalGrid,
    gcWeightFunction::GCWeightFunction,
    ConvolvedDensity::AbstractConvolvedDensity)
    for (zidx, zvalue) in enumerate(cosmologicalGrid.ZArray)
        gcWeightFunction.BiasArray[:, zidx] .=
        ComputeBias(zvalue, gcWeightFunction.BiasKind, ConvolvedDensity)
    end
end

function CreateBias(BiasDict::Dict)
    if BiasDict["model"] == "EuclidBias"
        Bias = CreateEuclidBias((BiasDict))
    else
        error("No Bias!")
    end
    return Bias
end

function CreateEuclidBias(BiasDict::Dict)
    Bias = CosmoCentral.EuclidBias()
    Bias.A = BiasDict["A"]
    Bias.B = BiasDict["B"]
    Bias.C = BiasDict["C"]
    Bias.D = BiasDict["D"]
    return Bias
end