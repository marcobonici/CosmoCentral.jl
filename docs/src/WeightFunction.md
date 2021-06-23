# Weight Function

## Galaxy Clustering
The expression of the Galaxy Clustering Weight Function is given by:

```math
W_{i}^{\mathrm{G}}(z)=b_{i}(z) n_{i}(z) \frac{H(z)}{c}
```

```@docs
CosmoCentral.GCWeightFunction
CosmoCentral.ComputeWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::CosmoCentral.AsbtractConvolvedDensity,
    AnalitycalDensity::CosmoCentral.AnalitycalDensityStruct,
    InstrumentResponse::CosmoCentral.InstrumentResponse,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmologyStruct,
    GCWeightFunction::CosmoCentral.GCWeightFunctionStruct)
CosmoCentral.ComputeWeightFunctionOverGrid(
    GCWeightFunction::CosmoCentral.GCWeightFunctionStruct,
    ConvolvedDensity::CosmoCentral.AsbtractConvolvedDensity,
    CosmologicalGrid::CosmoCentral.CosmologicalGrid,
    BackgroundQuantities::CosmoCentral.BackgroundQuantities,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmologyStruct)
```

# Weak Lensing
The expression of the Weak Lensing Weight Function is given by:

```math
W_{i}^{\gamma}(z)=\frac{3}{2}\left(\frac{H_{0}}{c}\right)^{2} \Omega_{\mathrm{m}, 0}(1+z) r(z) \widetilde{W}_{i}(z)
```

where ``\widetilde{W}_{i}(z)`` is the Lensing Efficiency, whose expression is given by

```math
\widetilde{W}_{i}(z)=\int_{z}^{z_{\max }} \mathrm{d} z^{\prime} n_{i}\left(z^{\prime}\right)\left[1-\frac{\tilde{r}(z)}{\tilde{r}\left(z^{\prime}\right)}\right]
```

```@docs
CosmoCentral.WLWeightFunction
CosmoCentral.ComputeWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::CosmoCentral.AsbtractConvolvedDensity,
    AnalitycalDensity::CosmoCentral.AnalitycalDensity,
    InstrumentResponse::CosmoCentral.InstrumentResponse,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmology,
    CosmologicalGrid::CosmoCentral.CosmologicalGrid,
    WLWeightFunction::CosmoCentral.WLWeightFunction)
CosmoCentral.ComputeWeightFunctionOverGrid(
    WLWeightFunction::CosmoCentral.WLWeightFunction,
    ConvolvedDensity::CosmoCentral.AsbtractConvolvedDensity,
    CosmologicalGrid::CosmoCentral.CosmologicalGrid,
    BackgroundQuantities::CosmoCentral.BackgroundQuantities,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmology)
CosmoCentral.ComputeLensingEfficiency(z::Float64, i::Int64,
    ConvolvedDensity::CosmoCentral.AsbtractConvolvedDensity,
    AnalitycalDensity::CosmoCentral.AnalitycalDensity,
    InstrumentResponse::CosmoCentral.InstrumentResponse,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmology,
    CosmologicalGrid::CosmoCentral.CosmologicalGrid,
    WLWeightFunction::CosmoCentral.WLWeightFunction)
CosmoCentral.ComputeLensingEfficiencyGrid!
CosmoCentral.ComputeLensingEfficiencyGrid!
```
