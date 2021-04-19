"""
    ComputeBias(z::Float64, ::PiecewiseBiasStruct,
    ConvolvedDensity::AsbtractConvolvedDensity)

This function evaluate the piecewise bias, given by ``\\sqrt{1+z}``, for a given
redshift.

"""
function ComputeBias(z::Float64, ::PiecewiseBiasStruct,
    ConvolvedDensity::AsbtractConvolvedDensity)
    idx = BinSearch(z, ConvolvedDensity.ZBinArray)
    bias = sqrt(1+(ConvolvedDensity.ZBinArray[idx]+
    ConvolvedDensity.ZBinArray[idx+1])/2)
    return bias
end

"""
    ComputeBiasOverGrid(CosmologicalGrid::CosmologicalGrid,
    GCWeightFunction::GCWeightFunctionStruct, Bias::AbstractBias,
    ConvolvedDensity::AsbtractConvolvedDensity)

This function evaluate the bias, over the cosmological redshift grid.
"""
function ComputeBiasOverGrid(CosmologicalGrid::CosmologicalGrid,
    GCWeightFunction::GCWeightFunctionStruct,
    ConvolvedDensity::AsbtractConvolvedDensity)
    for (zidx, zvalue) in enumerate(CosmologicalGrid.ZArray)
        GCWeightFunction.BiasArray[:, zidx] .=
        ComputeBias(zvalue, GCWeightFunction.BiasKind, ConvolvedDensity)
    end
end
