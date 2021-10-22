```@setup tutorial
using Plots; gr()
Plots.reset_defaults()
using CosmoCentral
using LaTeXStrings

w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()
CosmologicalGrid  = CosmoCentral.CosmologicalGrid(
ZArray=Array(LinRange(0.001, 4.0, 500)))

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
```@example tutorial

#instantiate background quantities and compute them
BackgroundQuantities = CosmoCentral.BackgroundQuantities(HZArray=
zeros(length(CosmologicalGrid.ZArray)),
χZArray=zeros(length(CosmologicalGrid.ZArray)))
CosmoCentral.ComputeBackgroundQuantitiesGrid!(CosmologicalGrid,
BackgroundQuantities, w0waCDMCosmology)

pH = plot(CosmologicalGrid.KArray, PowerSpectrum.PowerSpectrumNonlinArray[:,1],
ylabel = L"E(z)")
```