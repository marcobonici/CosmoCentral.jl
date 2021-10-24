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


CosmologicalGrid = CosmoCentral.CosmologicalGrid(ZArray = LinRange(0.001, 2.5, 500))
ConvolvedDensity = CosmoCentral.ConvolvedDensity(DensityGridArray =
        ones(10, length(CosmologicalGrid.ZArray)))
```

# Background Universe

Some of the most basics quantities in cosmology are the Hubble factor ``H(z)``
and the comoving distance ``\chi(z)``. In this section are presented the functions
which evaluate them.

## Hubble factor
```@docs
CosmoCentral.ComputeAdimensionalHubbleFactor
CosmoCentral.ComputeHubbleFactor
```

## Distances
```@docs
CosmoCentral.Computeχ
```

## BackgroundQuantities
In order to store them in a convenient way, we use a structure,
[`CosmoCentral.BackgroundQuantities`](@ref), which contains all the background quantities
```@docs
CosmoCentral.BackgroundQuantities
```

## Utils
A useful function, regarding the background quantities, is
[`CosmoCentral.ComputeBackgroundQuantitiesGrid!`](@ref). It computes the background
quantities over the redshift grid.
```@docs
CosmoCentral.ComputeBackgroundQuantitiesGrid!
```
Here we show how to calculate ``H(z)`` and ``\chi(z)`` over the redshift grid, then we plot 
them
```@example tutorial

#instantiate background quantities and compute them
BackgroundQuantities = CosmoCentral.BackgroundQuantities(HZArray=
zeros(length(CosmologicalGrid.ZArray)),
χZArray=zeros(length(CosmologicalGrid.ZArray)))
CosmoCentral.ComputeBackgroundQuantitiesGrid!(CosmologicalGrid,
BackgroundQuantities, w0waCDMCosmology)

pH = plot(CosmologicalGrid.ZArray, BackgroundQuantities.HZArray./ w0waCDMCosmology.H0,
ylabel = L"E(z)")
pχ = plot(CosmologicalGrid.ZArray, BackgroundQuantities.χZArray,
ylabel = L"\chi(z)\, \left[\mathrm{Mpc}\right]", xlabel= L"z")
plot(pH, pχ, layout = (2, 1), legend = false)
```
We also want to show that this code is fast. In this particular example, we are going to compute
``H(z)`` and ``\chi(z)`` on a redshift grid of 500 points
```@example tutorial
BackgroundQuantities = CosmoCentral.BackgroundQuantities(HZArray=
zeros(length(CosmologicalGrid.ZArray)),
χZArray=zeros(length(CosmologicalGrid.ZArray)))
@benchmark CosmoCentral.ComputeBackgroundQuantitiesGrid!(CosmologicalGrid,
BackgroundQuantities, w0waCDMCosmology)
```