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
        tickfontsize=12, legendfontsize=8, size = (550, 400), labelfontsize = 13,
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

CosmoCentral.ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
CosmoCentral.InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid, BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())


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
CosmoCentral.ComputeBiasGrid!(CosmologicalGrid, GCWeightFunction, ConvolvedDensity)
CosmoCentral.ComputeWeightFunctionGrid!(GCWeightFunction, ConvolvedDensity, CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)

WLWeightFunction = CosmoCentral.WLWeightFunction(WeightFunctionArray = zeros(length(ConvolvedDensity.DensityNormalizationArray), length(CosmologicalGrid.ZArray)), LensingEfficiencyArray = zeros(length(ConvolvedDensity.DensityNormalizationArray), length(CosmologicalGrid.ZArray)))
CosmoCentral.ComputeLensingEfficiencyGrid!(
    WLWeightFunction, ConvolvedDensity,
    CosmologicalGrid,
    BackgroundQuantities,
    w0waCDMCosmology, CosmoCentral.CustomLensingEfficiency())
CosmoCentral.ComputeWeightFunctionGrid!(WLWeightFunction, ConvolvedDensity, CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)

CosmoCentral.ComputeIntrinsicAlignmentGrid!(CosmologicalGrid, WLWeightFunction, ConvolvedDensity, BackgroundQuantities, w0waCDMCosmology)
CosmoCentral.ComputeWeightFunctionGrid!(WLWeightFunction, ConvolvedDensity, CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)

CℓGG = CosmoCentral.Cℓ(CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters), length(WLWeightFunction.WeightFunctionArray[:, 1]), length(WLWeightFunction.WeightFunctionArray[:, 1])))
CosmoCentral.ComputeCℓ!(CℓGG, GCWeightFunction, GCWeightFunction, BackgroundQuantities, w0waCDMCosmology, CosmologicalGrid, PowerSpectrum, CosmoCentral.CustomSimpson())

CℓGL = CosmoCentral.Cℓ(CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters), length(WLWeightFunction.WeightFunctionArray[:, 1]), length(WLWeightFunction.WeightFunctionArray[:, 1])))
CosmoCentral.ComputeCℓ!(CℓGL, GCWeightFunction, WLWeightFunction, BackgroundQuantities, w0waCDMCosmology, CosmologicalGrid, PowerSpectrum, CosmoCentral.CustomSimpson())
```



# Angular Coefficients
Central quantities are the Angular Coefficients ``C_\ell``. Actually we
implement only the Limber approximation to evaluate the ``C_\ell``,
[according to](https://arxiv.org/abs/1910.09273):

```math
C_{i j}^{AB}(\ell)=\frac{c}{H_0} \int_{z_{\min }}^{z_{\max }} \mathrm{d} z \frac{W_{i}^{A}(z) W_{j}^{B}(z)}{E(z) r^{2}(z)} P_{\delta \delta}\left(\frac{\ell+1 / 2}{r(z)}, z\right)
```

```@docs
CosmoCentral.Cℓ
CosmoCentral.∂Cℓ
CosmoCentral.ComputeCℓ!(Cℓ::CosmoCentral.AbstractCℓ,
WeightFunctionA::CosmoCentral.AbstractWeightFunction,
WeightFunctionB::CosmoCentral.AbstractWeightFunction,
BackgroundQuantities::CosmoCentral.BackgroundQuantities, ::CosmoCentral.w0waCDMCosmology,
CosmologicalGrid::CosmoCentral.AbstractCosmologicalGrid,
PowerSpectrum::CosmoCentral.AbstractPowerSpectrum,
::CosmoCentral.NumericalIntegrationSimpson)
CosmoCentral.ComputeCℓ!(Cℓ::CosmoCentral.AbstractCℓ,
WeightFunctionA::CosmoCentral.AbstractWeightFunction,
WeightFunctionB::CosmoCentral.AbstractWeightFunction,
BackgroundQuantities::CosmoCentral.BackgroundQuantities, ::CosmoCentral.w0waCDMCosmology,
CosmologicalGrid::CosmoCentral.AbstractCosmologicalGrid,
PowerSpectrum::CosmoCentral.AbstractPowerSpectrum, ::CosmoCentral.CustomSimpson)
```
Here we show how to compute and plot the  ``C_\ell``'s for Weak Lensing.
```@example tutorial
CℓLL = CosmoCentral.Cℓ(CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters),
length(WLWeightFunction.WeightFunctionArray[:, 1]),
length(WLWeightFunction.WeightFunctionArray[:, 1])))
CosmoCentral.ComputeCℓ!(CℓLL, WLWeightFunction, WLWeightFunction, BackgroundQuantities,
w0waCDMCosmology, CosmologicalGrid, PowerSpectrum, CosmoCentral.CustomSimpson())

x = CosmologicalGrid.ℓBinCenters
p = Plots.plot(xlabel=L"\ell", ylabel=L"\ell(\ell+1)C_{ii}^{LL}",
    title="Weak Lensing", legend=:bottomright)
for i in 1:10
    y =
    CℓLL.CℓArray[:,i,i] .* CosmologicalGrid.ℓBinCenters .*(CosmologicalGrid.ℓBinCenters .+1)
Plots.plot!(p, x, y, labels=(L"i=%$i"),  linewidth=3, xaxis=:log, yaxis=:log)
end
p
```
Following the same procedure it is possible to evaluate the Galaxy Clustering ``C_\ell``'s...
```@example tutorial
x = CosmologicalGrid.ℓBinCenters
p = Plots.plot(xlabel=L"\ell", ylabel=L"\ell(\ell+1)C_{ii}^{GG}",
    title="Galaxy Clustering", legend=:bottomright)
for i in 1:10
    y =
    CℓGG.CℓArray[:,i,i] .* CosmologicalGrid.ℓBinCenters .* (CosmologicalGrid.ℓBinCenters .+1)
Plots.plot!(p, x, y, labels=(L"i=%$i"),  linewidth=3, xaxis=:log, yaxis=:log)
end
p
```
...and the ``C_\ell``'s of the Cross-Correlation. Due to the contributions of the Intrinsic
Alignment, some of the ``C_\ell``'s are negative. Since we use the loglog
scale, we plot the absolute value of the angular coefficients; the negative ``C_\ell``'s are
indicated with dashed lines
```@example tutorial
x = CosmologicalGrid.ℓBinCenters
p = Plots.plot(xlabel=L"\ell", ylabel=L"\ell(\ell+1)|C_{ii}^{GL}|",
    title= "Galaxy-Galaxy Lensing", legend=:bottomright)
for i in 1:10
    y = CℓGL.CℓArray[:,i,i] .* CosmologicalGrid.ℓBinCenters .*
    (CosmologicalGrid.ℓBinCenters .+1)
    if any(x->x<=0, y)
        Plots.plot!(p, x, -y, labels=(L"i=%$i"),  linewidth=3, xaxis=:log, yaxis=:log,
        linestyle = :dash)
    else
        Plots.plot!(p, x, y, labels=(L"i=%$i"),  linewidth=3, xaxis=:log, yaxis=:log)
    end
end
p
```
Finally, the computation of the ``C_\ell``'s is quite fast. In particular, here we benchmark
their evaluation with 100 different ``\ell``'s value. 
```@example tutorial
@benchmark CosmoCentral.ComputeCℓ!(CℓGL, GCWeightFunction, WLWeightFunction, BackgroundQuantities, w0waCDMCosmology, CosmologicalGrid, PowerSpectrum, CosmoCentral.CustomSimpson())
```