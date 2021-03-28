function EvaluateTransferFunction(CosmologicalGrid::CosmologicalGridStruct,
    BackgroundQuantities::BackgroundQuantitiesStruct,
    κTransferFunction::κTransferFunctionStruct,
    PowerSpectrum::PowerSpectrumStruct)
    FFTLog = FFTLogStruct(XArray = CosmologicalGrid.ZArray, FXArray =
    κTransferFunction.WLWeightFunction.WeightFunctionArray .*
    PowerSpectrum.GrowthFactor)
    Kl , κTransferFunction.TransferFunctionArray =
    CosmoCentral.EvaluateFFTLog(FFTLog, CosmologicalGrid.MultipolesArray)
    κTransferFunction.TransferFunctionArray
    κTransferFunction.TransferFunctionArray =
    κTransferFunction.TransferFunctionArray ./ Kl ./ Kl .*
    CosmologicalGrid.MultipolesArray .* (CosmologicalGrid.MultipolesArray .+ 1)
end
