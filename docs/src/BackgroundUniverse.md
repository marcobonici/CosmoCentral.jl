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

## Utils
```@docs
CosmoCentral.ComputeBackgroundQuantitiesGrid!
```
Here we show how to calculate ``H(z)`` and ``\chi(z)`` over the redshift grid, then we plot 
them
```@example tutorial
using Plots
using CosmoCentral
w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()
CosmologicalGrid  = CosmoCentral.CosmologicalGrid(
ZArray=Array(LinRange(0.001, 4.0, 500)))
BackgroundQuantities = CosmoCentral.BackgroundQuantities(HZArray=
zeros(length(CosmologicalGrid.ZArray)),
χZArray=zeros(length(CosmologicalGrid.ZArray)))
CosmoCentral.ComputeBackgroundQuantitiesGrid!(CosmologicalGrid,
BackgroundQuantities, w0waCDMCosmology)
pH = plot(CosmologicalGrid.ZArray, BackgroundQuantities.HZArray./ w0waCDMCosmology.H0, ylabel = "E(z)", xlabel="z")
pχ = plot(CosmologicalGrid.ZArray, BackgroundQuantities.χZArray, ylabel = "χ(z) (Mpc)", xlabel="z")
plot(pH, pχ, layout = (2, 1), legend = false)
```
