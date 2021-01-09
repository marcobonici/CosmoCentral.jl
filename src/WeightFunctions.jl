"""
    ComputeGalaxyClusteringWeightFunction(z::Float64, i::Int64,
        ConvolvedDensity::ConvolvedDensity,
        w0waCDMCosmology::w0waCDMCosmology)

This function returns the source density for a given redshift ``z``.
"""
function ComputeWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::ConvolvedDensity, w0waCDMCosmology::w0waCDMCosmology,
    Bias::Bias, GCWeightFunction::GCWeightFunction)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    return ComputeConvolvedDensityFunction(z, i, ConvolvedDensity)*
    ComputeBias(z, Bias)
    ComputeHubbleFactor(z, w0waCDMCosmology) / c_0
end


"""
    ComputeGalaxyClusteringWeightFunction(CosmologicalGrid::CosmologicalGrid,
        ConvolvedDensity::ConvolvedDensity,
        w0waCDMCosmology::w0waCDMCosmology)

This function returns the source density for a given redshift ``z``.
"""
function ComputeWeightFunctionOverGrid(
    GCWeightFunction::GCWeightFunction, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, ConvolvedDensity::ConvolvedDensity,
    Bias::Bias, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    NormalizeConvolvedDensityStruct(ConvolvedDensity, AnalitycalDensity,
    InstrumentResponse)
    ComputeConvolvedDensityFunctionGrid(CosmologicalGrid, ConvolvedDensity,
    AnalitycalDensity, InstrumentResponse)
    ComputeBackgroundQuantitiesOverGrid(CosmologicalGrid, BackgroundQuantities,
    w0waCDMCosmology)
    ComputeBiasOverGrid(CosmologicalGrid, Bias, ConvolvedDensity)
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            GCWeightFunction.WeightFunctionArray[idx_ZBinArray, idx_ZArray] =
            Bias.BiasArray[idx_ZBinArray, idx_ZArray]*
            ConvolvedDensity.DensityGridArray[idx_ZBinArray,
            idx_ZArray] *
            BackgroundQuantities.HZArray[idx_ZArray] / c_0
        end
    end
end
