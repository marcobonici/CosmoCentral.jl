```@setup tutorial
using Plots; gr()
Plots.reset_defaults()
using PlotThemes
using CosmoCentral
using LaTeXStrings
using BenchmarkTools
default(palette = palette(:tab10))

w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()


plot_font = "Computer Modern"
Plots.default(titlefont = (16, plot_font), fontfamily=plot_font,
        linewidth=2, framestyle=:box, fg_legend =:black, label=nothing, grid=false,
        tickfontsize=12, legendfontsize=12, size = (550, 400), labelfontsize = 13,
        dpi = 200)

MultipolesArrayTemp = CosmoCentral.LogSpaced(10.,5000., 101)
MultipolesArray = zeros(length(MultipolesArrayTemp)-1)
#MultipolesWidths = vcat(CosmoCentral.Difference(MultipolesArrayTemp), ones(2000))
MultipolesWidths = CosmoCentral.Difference(MultipolesArrayTemp)
for i in 1:length(MultipolesWidths)
    MultipolesArray[i] = (MultipolesArrayTemp[i+1]+MultipolesArrayTemp[i])/2
end

path = joinpath(pwd(),"..","..","test","p_mm")
PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
CosmoCentral.ReadPowerSpectrumBackground(path, MultipolesArray, MultipolesWidths)
CosmoCentral.ExtractGrowthFactor!(BackgroundQuantities, PowerSpectrum)


AnalitycalDensity = CosmoCentral.AnalitycalDensity()
CosmoCentral.NormalizeAnalitycalDensity!(AnalitycalDensity)

#instantiate the instrument response and compute the convolved density
InstrumentResponse = CosmoCentral.InstrumentResponse()
ConvolvedDensity = CosmoCentral.ConvolvedDensity(DensityGridArray = ones(10,
length(CosmologicalGrid.ZArray)))
CosmoCentral.NormalizeConvolvedDensity!(ConvolvedDensity, AnalitycalDensity,
InstrumentResponse, CosmologicalGrid)
CosmoCentral.ComputeConvolvedDensityGrid!(CosmologicalGrid, ConvolvedDensity,
AnalitycalDensity, InstrumentResponse)
PiecewiseBias = CosmoCentral.PiecewiseBias()
GCWeightFunction = CosmoCentral.GCWeightFunction(WeightFunctionArray = zeros(length(ConvolvedDensity.DensityNormalizationArray), length(CosmologicalGrid.ZArray)))
WLWeightFunction = CosmoCentral.WLWeightFunction(WeightFunctionArray = zeros(length(ConvolvedDensity.DensityNormalizationArray), length(CosmologicalGrid.ZArray)), LensingEfficiencyArray = zeros(length(ConvolvedDensity.DensityNormalizationArray), length(CosmologicalGrid.ZArray)))
```

# Weight Function

## Galaxy Clustering
The expression of the Galaxy Clustering Weight Function is given by:

```math
W_{i}^{\mathrm{G}}(z)=b_{i}(z) n_{i}(z) \frac{H(z)}{c}
```

```@docs
CosmoCentral.GCWeightFunction
CosmoCentral.ComputeWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::CosmoCentral.AbstractConvolvedDensity,
    AnalitycalDensity::CosmoCentral.AnalitycalDensity,
    InstrumentResponse::CosmoCentral.InstrumentResponse,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmology,
    GCWeightFunction::CosmoCentral.GCWeightFunction)
CosmoCentral.ComputeWeightFunctionGrid!(
    GCWeightFunction::CosmoCentral.GCWeightFunction,
    ConvolvedDensity::CosmoCentral.AbstractConvolvedDensity,
    CosmologicalGrid::CosmoCentral.CosmologicalGrid,
    BackgroundQuantities::CosmoCentral.BackgroundQuantities,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmology)
```
For instance, here we plot the Galaxy Clustering weight function with a piecewise bias
```@example tutorial
PiecewiseBias = CosmoCentral.PiecewiseBias()
GCWeightFunction = CosmoCentral.GCWeightFunction(WeightFunctionArray =
zeros(length(ConvolvedDensity.DensityNormalizationArray), length(CosmologicalGrid.ZArray)))

CosmoCentral.ComputeBiasGrid!(CosmologicalGrid, GCWeightFunction, ConvolvedDensity)
CosmoCentral.ComputeWeightFunctionGrid!(GCWeightFunction, ConvolvedDensity,
CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
p = Plots.plot(xlabel=L"z", ylabel=L"W_i^g(z)\,\left[\mathrm{Mpc}^{-1}\right]")
for i in 1:10
Plots.plot!(p, CosmologicalGrid.ZArray, GCWeightFunction.WeightFunctionArray[i,:],
    labels=(L"i=%$i"),  linewidth=3)
end
p
```
Computing ``b(z)`` and ``W_i^g(z)`` is quite fast
```@example tutorial
@benchmark CosmoCentral.ComputeBiasGrid!(CosmologicalGrid, GCWeightFunction, ConvolvedDensity)
```

```@example tutorial
@benchmark CosmoCentral.ComputeWeightFunctionGrid!(GCWeightFunction, ConvolvedDensity,
CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
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
    ConvolvedDensity::CosmoCentral.AbstractConvolvedDensity,
    AnalitycalDensity::CosmoCentral.AnalitycalDensity,
    InstrumentResponse::CosmoCentral.InstrumentResponse,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmology,
    CosmologicalGrid::CosmoCentral.CosmologicalGrid,
    wlWeightFunction::CosmoCentral.WLWeightFunction)
CosmoCentral.ComputeWeightFunctionGrid!(
    wlWeightFunction::CosmoCentral.WLWeightFunction,
    ConvolvedDensity::CosmoCentral.AbstractConvolvedDensity,
    CosmologicalGrid::CosmoCentral.CosmologicalGrid,
    BackgroundQuantities::CosmoCentral.BackgroundQuantities,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmology)
CosmoCentral.ComputeLensingEfficiency(z::Float64, i::Int64,
    ConvolvedDensity::CosmoCentral.AbstractConvolvedDensity,
    AnalitycalDensity::CosmoCentral.AnalitycalDensity,
    InstrumentResponse::CosmoCentral.InstrumentResponse,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmology,
    CosmologicalGrid::CosmoCentral.CosmologicalGrid,
    wlWeightFunction::CosmoCentral.WLWeightFunction)
CosmoCentral.ComputeLensingEfficiencyGrid!(wlWeightFunction::CosmoCentral.WLWeightFunction,
    AnalitycalDensity::CosmoCentral.AnalitycalDensity,
    InstrumentResponse::CosmoCentral.InstrumentResponse,
    ConvolvedDensity::CosmoCentral.AbstractConvolvedDensity,
    CosmologicalGrid::CosmoCentral.CosmologicalGrid,
    BackgroundQuantities::CosmoCentral.BackgroundQuantities,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmology,
    ::CosmoCentral.StandardLensingEfficiency)
CosmoCentral.ComputeLensingEfficiencyGrid!(wlWeightFunction::CosmoCentral.WLWeightFunction,
    ConvolvedDensity::CosmoCentral.AbstractConvolvedDensity,
    CosmologicalGrid::CosmoCentral.CosmologicalGrid,
    BackgroundQuantities::CosmoCentral.BackgroundQuantities,
    w0waCDMCosmology::CosmoCentral.w0waCDMCosmology, ::CosmoCentral.CustomLensingEfficiency)
```
Here we plot the Weak Lensing weight function. In particular, the solid lines are pure shear,
while the dashed lines includes the Intrinsic Alignment contribution.
```@example tutorial
WLWeightFunction = CosmoCentral.WLWeightFunction(WeightFunctionArray = zeros(length(ConvolvedDensity.DensityNormalizationArray), length(CosmologicalGrid.ZArray)), LensingEfficiencyArray = zeros(length(ConvolvedDensity.DensityNormalizationArray), length(CosmologicalGrid.ZArray)))
CosmoCentral.ComputeLensingEfficiencyGrid!(
    WLWeightFunction, ConvolvedDensity,
    CosmologicalGrid,
    BackgroundQuantities,
    w0waCDMCosmology, CosmoCentral.CustomLensingEfficiency())
CosmoCentral.ComputeWeightFunctionGrid!(WLWeightFunction, ConvolvedDensity, CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
p = Plots.plot(xlabel=L"z", ylabel=L"W_i^\gamma(z)\,\left[\mathrm{Mpc}^{-1}\right]")
for i in 1:10
Plots.plot!(p, CosmologicalGrid.ZArray, WLWeightFunction.WeightFunctionArray[i,:],
    labels=(L"i=%$i"),  linewidth=3)
end

CosmoCentral.ComputeIntrinsicAlignmentGrid!(CosmologicalGrid, WLWeightFunction, ConvolvedDensity, BackgroundQuantities, w0waCDMCosmology)
CosmoCentral.ComputeWeightFunctionGrid!(WLWeightFunction, ConvolvedDensity, CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
for i in 1:10
Plots.plot!(p, CosmologicalGrid.ZArray, WLWeightFunction.WeightFunctionArray[i,:],
linewidth=3, linestyle = :dash)
end
p
```
Computations related to Weak Lensing are a bit slower, due to the nested integrals inside
``W_i^\gamma(z)``
```@example tutorial
@benchmark CosmoCentral.ComputeLensingEfficiencyGrid!(WLWeightFunction, ConvolvedDensity,
CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology,
CosmoCentral.CustomLensingEfficiency())
```