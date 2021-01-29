"""
    ComputeBias(z::Float64, PiecewiseBias::PiecewiseBias,
    ConvolvedDensity::ConvolvedDensity)

This function evaluate the piecewise bias, given by ``\\sqrt{1+z}``, for a given
redshift.

"""
function ComputeBias(z::Float64, PiecewiseBias::PiecewiseBiasStruct,
    ConvolvedDensity::AsbtractConvolvedDensity)
    idx = BinSearch(z, ConvolvedDensity.ZBinArray)
    bias = sqrt(1+(ConvolvedDensity.ZBinArray[idx]+
    ConvolvedDensity.ZBinArray[idx+1])/2)
    return bias
end

"""
    ComputeBiasOverGrid(CosmologicalGrid::CosmologicalGrid,
    PiecewiseBias::PiecewiseBiasStruct, ConvolvedDensity::ConvolvedDensity)

This function evaluate the piecewise bias, given by ``\\sqrt{1+z}``, over the
cosmological redshift grid.

"""
function ComputeBiasOverGrid(CosmologicalGrid::CosmologicalGrid,
    PiecewiseBias::PiecewiseBiasStruct, ConvolvedDensity::AsbtractConvolvedDensity)
    for (zidx, zvalue) in enumerate(CosmologicalGrid.ZArray)
        PiecewiseBias.BiasArray[:, zidx] .=
        ComputeBias(zvalue, PiecewiseBias, ConvolvedDensity)
    end
end
