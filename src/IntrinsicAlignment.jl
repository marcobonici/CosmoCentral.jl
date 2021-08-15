function ComputeIntrinsicAlignmentGrid!(CosmologicalGrid::CosmologicalGrid,
    LensingFunction::WLWeightFunction, ConvolvedDensity::AbstractConvolvedDensity,
    BackgroundQuantities::BackgroundQuantities, 
    Cosmology::AbstractCosmology)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    input_data = readdlm(joinpath(
        dirname(@__FILE__),"..","input_files","scaledmeanlum-E2Sa.txt"), Float64)
    z = input_data[:,1]
    lum = input_data[:,2]
    spl = Dierckx.Spline1D(z, lum, k = 2)
    lumgrid = spl.(CosmologicalGrid.ZArray)
    for zidx in 1:length(CosmologicalGrid.ZArray)
        for iidx in 1:length(ConvolvedDensity.ZBinArray)-1
            LensingFunction.IntrinsicAlignmentArray[iidx, zidx] =
            - BackgroundQuantities.HZArray[zidx] / c_0 *
            ConvolvedDensity.DensityGridArray[iidx, zidx] *
            LensingFunction.IntrinsicAlignmentModel.𝓐IA *
            LensingFunction.IntrinsicAlignmentModel.𝓒IA * Cosmology.ΩM *
            ( (1 + CosmologicalGrid.ZArray[zidx]) ^
            LensingFunction.IntrinsicAlignmentModel.ηIA ) *
            (lumgrid[zidx] ^ LensingFunction.IntrinsicAlignmentModel.βIA ) /
            BackgroundQuantities.DZArray[zidx]
        end
    end
end
