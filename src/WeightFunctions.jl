abstract type WeightFunction end

@kwdef mutable struct WeightFunctionStruct <: WeightFunction
    ConvolvedDensity::ConvolvedDensity = ConvolvedDensityStruct()
    PowerSpectrumGrid::PowerSpectrumGrid = PowerSpectrumGridStruct()
    w0waCDMCosmology::w0waCDMCosmology = w0waCDMStruct()
    Bias::Bias = PiecewiseBiasStruct()
    BackgroundQuantities::BackgroundQuantities = BackgroundQuantitiesStruct()
    WeightFunctionArray::AbstractArray{Float64, 2} =
    ones(length(ConvolvedDensity.zbinarray)-1, length(PowerSpectrumGrid.zgrid))
end

"""
    ComputeGalaxyClusteringWeightFunction(z::Float64, i::Int64,
        ConvolvedDensity::ConvolvedDensity,
        w0waCDMCosmology::w0waCDMCosmology)

This function returns the source density for a given redshift ``z``.
"""
function ComputeGalaxyClusteringWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::ConvolvedDensity, w0waCDMCosmology::w0waCDMCosmology,
    Bias::Bias)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    return ComputeConvolvedDensityFunction(z, i, ConvolvedDensity)*
    ComputeBias(z, Bias)
    ComputeHubbleFactor(z, w0waCDMCosmology) / c_0
end


"""
    ComputeGalaxyClusteringWeightFunction(PowerSpectrumGrid::PowerSpectrumGrid,
        ConvolvedDensity::ConvolvedDensity,
        w0waCDMCosmology::w0waCDMCosmology)

This function returns the source density for a given redshift ``z``.
"""
function ComputeGalaxyClusteringWeightFunctionOverGrid(
    WeightFunction::WeightFunction)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for idx_zbinarray in 1:length(WeightFunction.ConvolvedDensity.zbinarray)-1
        for idx_zgrid in 1:length(WeightFunction.PowerSpectrumGrid.zgrid)
            WeightFunction.WeightFunctionArray[idx_zbinarray, idx_zgrid] =
            WeightFunction.Bias.BiasArray[idx_zbinarray, idx_zgrid]*
            WeightFunction.ConvolvedDensity.densitygridarray[idx_zbinarray, idx_zgrid] *
            WeightFunction.BackgroundQuantities.Hzgrid[idx_zgrid] / c_0
        end
    end
end
