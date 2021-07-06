"""
    ComputeBias(z::Float64, ::PiecewiseBias,
    ConvolvedDensity::AbstractConvolvedDensity)

This function evaluate the piecewise bias, given by ``\\sqrt{1+z}``, for a given
redshift.

"""
function ComputeBias(z::Float64, Piecewise::PiecewiseBias,
    ConvolvedDensity::AbstractConvolvedDensity)
    idx = BinSearch(z, ConvolvedDensity.ZBinArray)
    bias = sqrt(1+(ConvolvedDensity.ZBinArray[idx]+
    ConvolvedDensity.ZBinArray[idx+1])/2) * Piecewise.BiasMultiplier[idx]
    return bias
end

function ComputeBias(z::Float64, Bias::EuclidBias,
    ConvolvedDensity::AbstractConvolvedDensity)
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
