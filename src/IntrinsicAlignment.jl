function ComputeIntrinsicAlignmentOverGrid(CosmologicalGrid::CosmologicalGrid,
    WLWeightFunction::WLWeightFunctionStruct,
    ConvolvedDensity::AbstractConvolvedDensity,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    input_data = readdlm("input_files/scaledmeanlum-E2Sa.txt", Float64)
    z = input_data[:,1]
    lum = input_data[:,2]
    spl = Dierckx.Spline1D(z, lum)
    for (zidx, zvalue) in enumerate(CosmologicalGrid.ZArray)
        for iidx in 1:length(ConvolvedDensity.ZBinArray)-1
            WLWeightFunction.IntrinsicAlignmentArray[iidx, zidx] =
            -BackgroundQuantities.HZArray[zidx]/c_0 *
            ConvolvedDensity.DensityGridArray[iidx, zidx] *
            WLWeightFunction.IntrinsicAlignmentModel.A *
            WLWeightFunction.IntrinsicAlignmentModel.C *
            w0waCDMCosmology.ΩM *
            (1 + zidx) ^ WLWeightFunction.IntrinsicAlignmentModel.η *
            spl(zvalue) ^ WLWeightFunction.IntrinsicAlignmentModel.β
        end
    end
end
