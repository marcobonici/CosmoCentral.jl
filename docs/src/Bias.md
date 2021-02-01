```@setup tutorial
using Plots CosmoCentral
pyplot()
import PyPlot
```

# Bias

One of the most important ingredients of the Galaxy Clustering Weight function
is the Galaxy Bias, that is, the statistical relation between the
distribution of galaxies and matter

```math
\delta_{g}(\boldsymbol{x})=b \delta(\boldsymbol{x})
```

In this page are listed the Bias Model we have
implemented.

## PiecewiseBias

The Piecewise Bias model is taken from the
[official Euclid forecast](https://arxiv.org/abs/1910.09273)

```math
b(z)= \sqrt{1+\bar{z}}
```

where ``\bar{z}`` is the redshift value in the center of the tomographic bin
where the redshift ``z`` lies.
The PiecewiseBias here is plotted.
```@example tutorial
CosmologicalGrid = CosmoCentral.CosmologicalGridStruct(
        ZArray = LinRange(0.001, 2.5, 500))
ConvolvedDensity = CosmoCentral.ConvolvedDensityStruct(DensityGridArray =
        ones(10, length(CosmologicalGrid.ZArray)))
GCWeightFunction = CosmoCentral.GCWeightFunctionStruct(WeightFunctionArray=
        zeros(length(ConvolvedDensity.DensityGridArray[1,:]),
        length(CosmologicalGrid.ZArray)))
CosmoCentral.ComputeBiasOverGrid(CosmologicalGrid, GCWeightFunction,
GCWeightFunction.BiasKind,
ConvolvedDensity)
x = CosmologicalGrid.ZArray; y = GCWeightFunction.BiasArray[1, :];
plot(x, y, label = "Bias", xlabel="z")
```

```@docs
CosmoCentral.ComputeBias
CosmoCentral.ComputeBiasOverGrid
```
