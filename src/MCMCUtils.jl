function EvaluateCℓMCMCStep(Cosmology::AbstractCosmology, Density::AbstractConvolvedDensity,
    Bias::AbstractBias, IA::AbstractIntrinsicAlignment, CosmoGrid::AbstractCosmologicalGrid)
    backgroundquantities = BackgroundQuantities(HZArray= zeros(length(CosmoGrid.ZArray)),
    rZArray=zeros(length(CosmoGrid.ZArray)))
    ComputeBackgroundQuantitiesGrid!(CosmoGrid, backgroundquantities, Cosmology)
    ClassyParams = Initializeclassy(Cosmology)
    Pmm = PowerSpectrum(PowerSpectrumLinArray =
    zeros(length(CosmoGrid.KArray), length(CosmoGrid.ZArray)),
    PowerSpectrumNonlinArray = zeros(length(CosmoGrid.KArray), length(CosmoGrid.ZArray)),
    InterpolatedPowerSpectrum = zeros(length(CosmoGrid.ℓBinCenters),
    length(CosmoGrid.ZArray)))
    EvaluatePowerSpectrum!(ClassyParams, CosmoGrid, Pmm)
    ComputeLimberArray!(CosmoGrid, backgroundquantities)
    InterpolatePowerSpectrumLimberGrid!(CosmoGrid, backgroundquantities, Pmm,
    BSplineCubic())
    ExtractGrowthFactor!(backgroundquantities, Pmm)
    GCW = GCWeightFunction(WeightFunctionArray =
    zeros(length(Density.DensityNormalizationArray), length(CosmoGrid.ZArray)),
    BiasKind = Bias)
    CosmoCentral.ComputeBiasGrid!(CosmoGrid, GCW, Density)
    CosmoCentral.ComputeWeightFunctionGrid!(GCW, Density, CosmoGrid,
    backgroundquantities, Cosmology)
    WLW = WLWeightFunction(WeightFunctionArray =
    zeros(length(Density.DensityNormalizationArray), length(CosmoGrid.ZArray)),
    LensingEfficiencyArray = zeros(length(Density.DensityNormalizationArray),
    length(CosmoGrid.ZArray)), IntrinsicAlignmentModel = IA)
    ComputeLensingEfficiencyGrid!(WLW, Density, CosmoGrid,
    backgroundquantities, Cosmology, CustomLensingEfficiency())
    ComputeIntrinsicAlignmentGrid!(CosmoGrid, WLW, Density,
    backgroundquantities, Cosmology, "../inputs/scaledmeanlum-E2Sa.txt")
    #TODO pay attention at this path!
    ComputeWeightFunctionGrid!(WLW, Density, CosmoGrid, backgroundquantities,
    Cosmology)
    CℓLL = Cℓ(CℓArray = zeros(length(CosmoGrid.ℓBinCenters),
    length(WLW.WeightFunctionArray[:, 1]),
    length(WLW.WeightFunctionArray[:, 1])))
    ComputeCℓ!(CℓLL, WLW, WLW, backgroundquantities, Cosmology,
    CosmoGrid, Pmm, CustomSimpson())
    return CℓLL
end