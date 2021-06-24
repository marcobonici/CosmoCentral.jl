function ComputeIntrinsicAlignmentGrid!(CosmologicalGrid::CosmologicalGrid,
    LensingFunction::WLWeightFunction, ConvolvedDensity::AbstractConvolvedDensity,
    BackgroundQuantities::BackgroundQuantities, PowerSpectrum::PowerSpectrum,
    w0waCDMCosmology::w0waCDMCosmology, Path::String)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    input_data = readdlm(Path, Float64)
    z = input_data[:,1]
    lum = input_data[:,2]
    spl = Dierckx.Spline1D(z, lum, k = 1)
    for (zidx, zvalue) in enumerate(CosmologicalGrid.ZArray)
        for iidx in 1:length(ConvolvedDensity.ZBinArray)-1
            LensingFunction.IntrinsicAlignmentArray[iidx, zidx] =
            - BackgroundQuantities.HZArray[zidx] / c_0 *
            ConvolvedDensity.DensityGridArray[iidx, zidx] *
            LensingFunction.IntrinsicAlignmentModel.A *
            LensingFunction.IntrinsicAlignmentModel.C * w0waCDMCosmology.ΩM *
            ( (1 + zvalue) ^ LensingFunction.IntrinsicAlignmentModel.η ) *
            (spl(zvalue) ^ LensingFunction.IntrinsicAlignmentModel.β ) /
            PowerSpectrum.GrowthFactor[zidx]
        end
    end
end
