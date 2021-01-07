function ComputeBias(z::Float64, PiecewiseBias::PiecewiseBias)
    idx = BinSearch(z, PiecewiseBias.ConvolvedDensity.ZBinArray)
    bias = sqrt(1+(PiecewiseBias.ConvolvedDensity.ZBinArray[idx]+
    PiecewiseBias.ConvolvedDensity.ZBinArray[idx+1])/2)
    return bias
end

function ComputeBiasOverGrid(PiecewiseBias::PiecewiseBias)
    for (zidx, zvalue) in enumerate(PiecewiseBias.ConvolvedDensity.ZBinArray)
        PiecewiseBias.BiasArray[:, zidx] .=
        ComputeBias(zvalue, PiecewiseBias)
    end
end
