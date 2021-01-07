function ComputeBias(z::Float64, PiecewiseBias::PiecewiseBias,
    ConvolvedDensity::ConvolvedDensity)
    idx = BinSearch(z, ConvolvedDensity.ZBinArray)
    bias = sqrt(1+(ConvolvedDensity.ZBinArray[idx]+
    ConvolvedDensity.ZBinArray[idx+1])/2)
    return bias
end

function ComputeBiasOverGrid(CosmologicalGrid::CosmologicalGrid,
    PiecewiseBias::PiecewiseBias, ConvolvedDensity::ConvolvedDensity)
    for (zidx, zvalue) in enumerate(CosmologicalGrid.ZArray)
        PiecewiseBias.BiasArray[:, zidx] .=
        ComputeBias(zvalue, PiecewiseBias, ConvolvedDensity)
    end
end
