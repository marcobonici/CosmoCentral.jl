abstract type Bias end
abstract type PiecewiseBias <: Bias end

@kwdef mutable struct PiecewiseBiasStruct <: PiecewiseBias
    ConvolvedDensity::ConvolvedDensity = ConvolvedDensityStruct()
    PowerSpectrumGrid::PowerSpectrumGrid = PowerSpectrumGridStruct()
    BiasArray::AbstractArray{Float64, 2} =
    ones(length(ConvolvedDensity.zbinarray)-1, length(PowerSpectrumGrid.zgrid))
end

function ComputeBias(z::Float64, PiecewiseBias::PiecewiseBias)
    idx = BinSearch(z, PiecewiseBias.ConvolvedDensity.zbinarray)
    bias = sqrt(1+(PiecewiseBias.ConvolvedDensity.zbinarray[idx]+
    PiecewiseBias.ConvolvedDensity.zbinarray[idx+1])/2)
    return bias
end

function ComputeBiasOverGrid(PiecewiseBias::PiecewiseBias)
    for (zidx, zvalue) in enumerate(PiecewiseBias.ConvolvedDensity.zbinarray)
        PiecewiseBias.BiasArray[:, zidx] .=
        ComputeBias(zvalue, PiecewiseBias)
    end
end
