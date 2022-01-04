function ComputeSourceFunctionOverGrid(
    LensingSourceFunction::LensingSourceFunction,
    ConvolvedDensity::AbstractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    Cosmology::AbstractCosmology,
    PowerSpectrum::PowerSpectrum)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for zbinidx in 1:length(ConvolvedDensity.ZBinArray)-1
        for zidx in 1:length(CosmologicalGrid.ZArray)
            LensingSourceFunction.SourceFunctionArray[zbinidx, zidx] =
            1.5 * (Cosmology.H0/c_0)^2 * Cosmology.ΩM *
            (1. + CosmologicalGrid.ZArray[zidx]) *
            LensingSourceFunction.LensingEfficiencyArray[zbinidx, zidx] *
            BackgroundQuantities.χZArray[zidx]
        end
    end
end
