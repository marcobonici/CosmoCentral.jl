function ComputeIntrinsicAlignmentGrid!(CosmologicalGrid::CosmologicalGrid,
    LensingFunction::WLWeightFunction, ConvolvedDensity::AbstractConvolvedDensity,
    BackgroundQuantities::BackgroundQuantities, 
    Cosmology::AbstractCosmology, Path::String)
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
            LensingFunction.IntrinsicAlignmentModel.𝓐IA *
            LensingFunction.IntrinsicAlignmentModel.𝓒IA * Cosmology.ΩM *
            ( (1 + zvalue) ^ LensingFunction.IntrinsicAlignmentModel.ηIA ) *
            (spl(zvalue) ^ LensingFunction.IntrinsicAlignmentModel.βIA ) /
            BackgroundQuantities.DZArray[zidx]
        end
    end
end
