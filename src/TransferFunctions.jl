function EvaluateTransferFunction(CosmologicalGrid::CosmologicalGridStruct,
    BackgroundQuantities::BackgroundQuantitiesStruct,
    ConvolvedDensity::ConvolvedDensityStruct,
    κTransferFunction::κTransferFunctionStruct,
    PowerSpectrum::PowerSpectrumStruct)
    c_0 = 2.99792458e5
    GrowthFactorZ = Dierckx.Spline1D(CosmologicalGrid.ZArray,
    PowerSpectrum.GrowthFactor)
    for iidx in 1:length(ConvolvedDensity.ZBinArray)-1
        FFTLog = FFTLogStruct(XArray = BackgroundQuantities.rZArray, FXArray =
        κTransferFunction.LensingSourceFunction.SourceFunctionArray[iidx, :] .*
        BackgroundQuantities.rZArray .* PowerSpectrum.GrowthFactor ./
        GrowthFactorZ((
        ConvolvedDensity.ZBinArray[iidx+1]-ConvolvedDensity.ZBinArray[iidx])/2))
        Kl, κTransferFunction.TransferFunctionArray[iidx, :, :] =
        CosmoCentral.EvaluateFFTLog(FFTLog, CosmologicalGrid.MultipolesArray)
        κTransferFunction.TransferFunctionArray[iidx, :, :] =
        κTransferFunction.TransferFunctionArray[iidx, :, :] ./
        Kl ./ Kl .* CosmologicalGrid.MultipolesArray .*
        (CosmologicalGrid.MultipolesArray .+ 1)
        CosmologicalGrid.KBeyondLimberArray = Kl
    end
end
