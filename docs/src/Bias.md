```@setup tutorial
using Plots, CosmoCentral
pyplot()
import PyPlot
```

# Bias

One of the most important ingredients of the Galaxy Clustering Weight function
is the Galaxy Bias. In this page are listed the Bias Model we have
implemented.

## PiecewiseBias

The Piecewise Bias is taken from the
[official Euclid forecast](https://arxiv.org/abs/1910.09273)
```math
b(z)= \sqrt{1+z}
```
and is evaluated in the centre of the tomographic bins.
The bias here is initialized and plotted.
```@example tutorial
CosmologicalGrid = CosmoCentral.CosmologicalGridStruct(
        ZArray = LinRange(0.001, 2.5, 300))
ConvolvedDensity = CosmoCentral.ConvolvedDensityStruct(DensityGridArray =
        ones(10, length(CosmologicalGrid.ZArray)))
Bias = CosmoCentral.PiecewiseBiasStruct(
        BiasArray = zeros(length(ConvolvedDensity.DensityNormalizationArray),
        length(CosmologicalGrid.ZArray)))
CosmoCentral.ComputeBiasOverGrid(CosmologicalGrid, Bias, ConvolvedDensity)
x = CosmologicalGrid.ZArray; y = Bias.BiasArray[1, :];
plot(x, y, label = "Bias", xlabel="z")
```

```@docs
CosmoCentral.ComputeBias
CosmoCentral.ComputeBiasOverGrid
```
