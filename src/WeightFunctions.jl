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
    ComputeGalaxyClusteringWeightFunction(CosmologicalGrid::CosmologicalGrid,
        ConvolvedDensity::ConvolvedDensity,
        w0waCDMCosmology::w0waCDMCosmology)

This function returns the source density for a given redshift ``z``.
"""
function ComputeGalaxyClusteringWeightFunctionOverGrid(
    WeightFunction::WeightFunction, ConvolvedDensity::ConvolvedDensity,
    Bias::Bias)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            WeightFunction.WeightFunctionArray[idx_ZBinArray, idx_ZArray] =
            Bias.BiasArray[idx_ZBinArray, idx_ZArray]*
            ConvolvedDensity.DensityGridArray[idx_ZBinArray,
            idx_ZArray] *
            WeightFunction.BackgroundQuantities.HZArray[idx_ZArray] / c_0
        end
    end
end
