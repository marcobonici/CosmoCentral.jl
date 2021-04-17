function ComputeSourceFunctionOverGrid(
    LensingSourceFunction::LensingSourceFunctionStruct,
    ConvolvedDensity::AsbtractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    PowerSpectrum::PowerSpectrumStruct)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            LensingSourceFunction.SourceFunctionArray[idx_ZBinArray, idx_ZArray] =
            1.5 * (w0waCDMCosmology.H0/c_0)^2 * w0waCDMCosmology.ΩM *
            (1. + CosmologicalGrid.ZArray[idx_ZArray]) *
            LensingSourceFunction.LensingEfficiencyArray[idx_ZBinArray, idx_ZArray] *
            PowerSpectrum.GrowthFactor[idx_ZArray] /
            PowerSpectrum.GrowthFactor[1]
        end
    end
end
