function EvaluateTransferFunction(CosmologicalGrid::CosmologicalGridStruct,
    BackgroundQuantities::BackgroundQuantitiesStruct,
    ConvolvedDensity::ConvolvedDensityStruct,
    κTransferFunction::κTransferFunctionStruct,
    PowerSpectrum::PowerSpectrumStruct)
    for iidx in 1:length(ConvolvedDensity.ZBinArray)-1
        FFTLog = FFTLogStruct(XArray = BackgroundQuantities.rZArray, FXArray =
        κTransferFunction.LensingSourceFunction.SourceFunctionArray[iidx, :] .*
        PowerSpectrum.GrowthFactor)
        Kl, κTransferFunction.TransferFunctionArray[iidx, :, :] =
        CosmoCentral.EvaluateFFTLog(FFTLog, CosmologicalGrid.MultipolesArray)
        κTransferFunction.TransferFunctionArray[iidx, :, :] =
        κTransferFunction.TransferFunctionArray[iidx, :, :] ./
        Kl ./ Kl .* CosmologicalGrid.MultipolesArray .*
        (CosmologicalGrid.MultipolesArray .+ 1)
        CosmologicalGrid.KBeyondLimberArray = Kl
    end
end
