function ComputeSourceFunctionOverGrid(
    LensingSourceFunction::LensingSourceFunctionStruct,
    ConvolvedDensity::AbstractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    PowerSpectrum::PowerSpectrumStruct)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for zbinidx in 1:length(ConvolvedDensity.ZBinArray)-1
        for zidx in 1:length(CosmologicalGrid.ZArray)
            LensingSourceFunction.SourceFunctionArray[zbinidx, zidx] =
            1.5 * (w0waCDMCosmology.H0/c_0)^2 * w0waCDMCosmology.ΩM *
            (1. + CosmologicalGrid.ZArray[zidx]) *
            LensingSourceFunction.LensingEfficiencyArray[zbinidx, zidx] *
            BackgroundQuantities.rZArray[zidx]
        end
    end
end
