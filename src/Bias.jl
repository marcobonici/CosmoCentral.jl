"""
    ComputeBias(z::Float64, ::PiecewiseBias,
    ConvolvedDensity::AbstractConvolvedDensity)

This function evaluate the piecewise bias, given by ``\\sqrt{1+z}``, for a given
redshift.

"""
function ComputeBias(z::T, Piecewise::PiecewiseBias,
    convdens::AbstractConvolvedDensity) where T
    idx = BinSearch(z, convdens.ZBinArray)
    bias = sqrt(1+(convdens.ZBinArray[idx]+
    convdens.ZBinArray[idx+1])/2) * Piecewise.BiasMultiplier[idx]
    return bias
end

function ComputeBias(z::T, Piecewise::PiecewiseBias,
    convdens::AbstractConvolvedDensity, input_array) where T
    idx = BinSearch(z, convdens.ZBinArray)
    bias = input_array[idx]
    return bias
end

function ComputeBias(z::T, Bias::EuclidBias,
    convdens::AbstractConvolvedDensity) where T
    bias = Bias.A+Bias.B/(1+exp((Bias.D-z)*Bias.C))
    return bias
end

"""
    ComputeBiasGrid!(cosmogrid::CosmologicalGrid,
    gcWeightFunction::GCWeightFunction, Bias::AbstractBias,
    ConvolvedDensity::AbstractConvolvedDensity)

This function evaluate the bias, over the cosmological redshift grid.
"""
function ComputeBiasGrid!(cosmogrid::CosmologicalGrid,GCW::GCWeightFunction,
    convdens::AbstractConvolvedDensity)
    for (zidx, zvalue) in enumerate(cosmogrid.ZArray)
        GCW.BiasArray[:, zidx] .=
        ComputeBias(zvalue, GCW.BiasKind, convdens)
    end
end
function ComputeBiasGrid!(cosmogrid::CosmologicalGrid,GCW::GCWeightFunction,
    convdens::AbstractConvolvedDensity, input_array)
    for (zidx, zvalue) in enumerate(cosmogrid.ZArray)
        GCW.BiasArray[:, zidx] .=
        ComputeBias(zvalue, GCW.BiasKind, convdens, input_array)
    end
end

function CreateBias(BiasDict::Dict)
    if BiasDict["model"] == "EuclidBias"
        Bias = CreateEuclidBias((BiasDict))
    elseif BiasDict["model"] == "PiecewiseBias"
        Bias = CreateISTFBias((BiasDict))
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

function CreateISTFBias(BiasDict::Dict)
    BiasMultiplier = [BiasDict["b1"], BiasDict["b2"], BiasDict["b3"], BiasDict["b4"], BiasDict["b5"], BiasDict["b6"], BiasDict["b7"], BiasDict["b8"], BiasDict["b9"], BiasDict["b10"]]
    Bias = CosmoCentral.PiecewiseBias()
    Bias.BiasMultiplier = BiasMultiplier
    return Bias
end

function CreateISTFBias(BiasDict::Dict)
    BiasMultiplier = [BiasDict["b1"], BiasDict["b2"], BiasDict["b3"], BiasDict["b4"], BiasDict["b5"], BiasDict["b6"], BiasDict["b7"], BiasDict["b8"], BiasDict["b9"], BiasDict["b10"]]
    Bias = CosmoCentral.PiecewiseBias()
    Bias.BiasMultiplier = BiasMultiplier
    return Bias
end
