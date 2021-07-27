function EvaluateCℓMCMCStep(Cosmology::AbstractCosmology, Density::AbstractConvolvedDensity,
    Bias::AbstractBias, IA::AbstractIntrinsicAlignment, CosmoGrid::AbstractCosmologicalGrid)
    backgroundquantities = BackgroundQuantities(HZArray= zeros(length(CosmoGrid.ZArray)),
    rZArray=zeros(length(CosmoGrid.ZArray)))
    ComputeBackgroundQuantitiesGrid!(CosmoGrid, backgroundquantities, Cosmology)
    ClassyParams = Initializeclassy(Cosmology)
    PowerSpectrum = PowerSpectrum(PowerSpectrumLinArray =
    zeros(length(CosmoGrid.KArray), length(CosmoGrid.ZArray)),
    PowerSpectrumNonlinArray = zeros(length(CosmoGrid.KArray), length(CosmoGrid.ZArray)),
    InterpolatedPowerSpectrum = zeros(length(CosmoGrid.ℓBinCenters), length(CosmoGrid.ZArray)))
    EvaluatePowerSpectrum!(ClassyParams, CosmoGrid, PowerSpectrum)
    ComputeLimberArray!(CosmoGrid, backgroundquantities)
    InterpolatePowerSpectrumLimberGrid!(CosmoGrid, backgroundquantities, PowerSpectrum,
    BSplineCubic())
    ExtractGrowthFactor!(backgroundquantities, PowerSpectrum)
    GCWeightFunction = GCWeightFunction(WeightFunctionArray =
    zeros(length(Density.DensityNormalizationArray), length(CosmoGrid.ZArray)), BiasKind = Bias)
    CosmoCentral.ComputeBiasGrid!(CosmoGrid, GCWeightFunction, Density)
    CosmoCentral.ComputeWeightFunctionGrid!(GCWeightFunction, Density, CosmoGrid,
    backgroundquantities, Cosmology)
    WLWeightFunction = WLWeightFunction(WeightFunctionArray =
    zeros(length(Density.DensityNormalizationArray), length(CosmoGrid.ZArray)),
    LensingEfficiencyArray = zeros(length(Density.DensityNormalizationArray),
    length(CosmoGrid.ZArray)))
    ComputeLensingEfficiencyGrid!(WLWeightFunction, Density, CosmoGrid,
    backgroundquantities, Cosmology, CustomLensingEfficiency())
    ComputeIntrinsicAlignmentGrid!(CosmoGrid, WLWeightFunction, Density,
    backgroundquantities, Cosmology, "../inputs/scaledmeanlum-E2Sa.txt")
    #TODO pay attention at this path!
    ComputeWeightFunctionGrid!(WLWeightFunction, Density, CosmoGrid, backgroundquantities,
    Cosmology)
    CℓLL = Cℓ(CℓArray = zeros(length(CosmoGrid.ℓBinCenters),
    length(WLWeightFunction.WeightFunctionArray[:, 1]),
    length(WLWeightFunction.WeightFunctionArray[:, 1])))
    ComputeCℓ!(CℓLL, WLWeightFunction, WLWeightFunction, backgroundquantities, Cosmology,
    CosmoGrid, PowerSpectrum, CustomSimpson())
    return CℓLL
end