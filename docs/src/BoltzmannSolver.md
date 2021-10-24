```@setup tutorial
using Plots; gr()
Plots.reset_defaults()
using CosmoCentral
using LaTeXStrings
using BenchmarkTools

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
CosmoCentral.ExtractGrowthFactor!(BackgroundQuantities, PowerSpectrum)

CosmoCentral.ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
```

# BoltzmannSolver

One of the core components of CosmoCentral are the Boltzmann solvers, codes
developed by the community that solve the coupled Einstein-Boltzmann equations
and used to evaluate the Matter Power Spectrum. Currently, the Boltzmann codes
implemented in CosmoCentral are:

- CLASS

We plan to include other Boltzmann solvers in the future.

## Power Spectrum

```@docs
CosmoCentral.PowerSpectrum
CosmoCentral.ComputeLimberArray!
CosmoCentral.InterpolatePowerSpectrumLimberGrid!
CosmoCentral.InterpolationMethod
```

## classy

classy is the Python wrapper for [CLASS](http://class-code.net/).

```@docs
CosmoCentral.classyParams
CosmoCentral.Initializeclassy
CosmoCentral.EvaluatePowerSpectrum!
```
For instance, here we show the linear and nonlinear ``P_{\delta\delta}(k,z)``, evaluated by
CLASS, for a the reference cosmology.
```@example tutorial
p = plot(CosmologicalGrid.KArray, PowerSpectrum.PowerSpectrumLinArray[:,1],
ylabel = L"P(k)\,\left[\mathrm{Mpc}^{3}\right]", xlabel = L"k\,\left[\mathrm{Mpc}^{-1}\right]",
xaxis=:log, yaxis=:log, label = L"\mathrm{Linear}", xlims = (1e-5,10), ylims = (1,1e5))
plot!(p, CosmologicalGrid.KArray, PowerSpectrum.PowerSpectrumNonlinArray[:,1],
xaxis=:log, yaxis=:log, label = L"\mathrm{Nonlinear}")
```
When performing forecasts, ``P_{\delta \delta}(k,z)`` can be valuated once and stored (we 
provide a precomputed set of spectra [here](https://zenodo.org/record/5270335)), so
their computational impact is reduced. However, when evaluating ``C_\ell``'s, there is a not
negligible impact of ``P_{ \delta\delta}(k,z)`` interpolation and evaluation on the Limber 
``k-``grid.
```@example tutorial
@benchmark CosmoCentral.InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
```