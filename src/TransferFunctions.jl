function EvaluateTransferFunction(CosmologicalGrid::CosmologicalGridStruct,
    BackgroundQuantities::BackgroundQuantitiesStruct,
    ConvolvedDensity::ConvolvedDensityStruct,
    κTransferFunction::κTransferFunctionStruct,
    PowerSpectrum::PowerSpectrumStruct)
    c_0 = 2.99792458e5
    for iidx in 1:length(ConvolvedDensity.ZBinArray)-1
        FFTLog = FFTLogStruct(XArray = CosmologicalGrid.ZArray, FXArray =
        κTransferFunction.WLWeightFunction.WeightFunctionArray[iidx, :] .*
        PowerSpectrum.GrowthFactor)
        #TODO : add division by GrowthFactor(z_i)
        Kl, κTransferFunction.TransferFunctionArray[iidx, :, :] =
        CosmoCentral.EvaluateFFTLog(FFTLog, CosmologicalGrid.MultipolesArray)
        κTransferFunction.TransferFunctionArray
        κTransferFunction.TransferFunctionArray =
        κTransferFunction.TransferFunctionArray ./ Kl ./ Kl ./ (c_0^2) .*
        CosmologicalGrid.MultipolesArray .* (CosmologicalGrid.MultipolesArray .+ 1)
        CosmologicalGrid.KBeyondLimberArray = Kl
    end
end
