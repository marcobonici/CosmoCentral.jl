function EvaluateTransferFunction(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    ConvolvedDensity::ConvolvedDensity,
    κTransferFunction::κTransferFunction,
    PowerSpectrum::PowerSpectrum)
    for iidx in 1:length(ConvolvedDensity.ZBinArray)-1
        fftlog = FFTLog(XArray = BackgroundQuantities.χZArray, FXArray =
        κTransferFunction.LensingSourceFunction.SourceFunctionArray[iidx, :] .*
        PowerSpectrum.GrowthFactor)
        Kl, κTransferFunction.TransferFunctionArray[iidx, :, :] =
        CosmoCentral.EvaluateFFTLog(fftlog, CosmologicalGrid.ℓBinCenters)
        κTransferFunction.TransferFunctionArray[iidx, :, :] =
        κTransferFunction.TransferFunctionArray[iidx, :, :] ./
        Kl ./ Kl .* CosmologicalGrid.ℓBinCenters .*
        (CosmologicalGrid.ℓBinCenters .+ 1)
        CosmologicalGrid.KBeyondLimberArray = Kl
    end
end
