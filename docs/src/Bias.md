```@setup tutorial
using Plots; gr()
Plots.reset_defaults()
using CosmoCentral
using LaTeXStrings

plot_font = "Computer Modern"
Plots.default(titlefont = (16, plot_font), fontfamily=plot_font,
        linewidth=2, framestyle=:box, fg_legend =:black, label=nothing, grid=false,
        tickfontsize=12, legendfontsize=12, size = (550, 400), labelfontsize = 13,
        dpi = 200)


CosmologicalGrid = CosmoCentral.CosmologicalGrid(ZArray = LinRange(0.001, 2.5, 500))
ConvolvedDensity = CosmoCentral.ConvolvedDensity(DensityGridArray =
        ones(10, length(CosmologicalGrid.ZArray)))
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

## Piecewise Bias

The Piecewise Bias model is taken from the
[official Euclid forecast](https://arxiv.org/abs/1910.09273)

```math
b(z)= \sqrt{1+\bar{z}}
```

where ``\bar{z}`` is the redshift value in the center of the tomographic bin
where the redshift ``z`` lies.
The PiecewiseBias is plotted here.
```@example tutorial
GCWeightFunction = CosmoCentral.GCWeightFunction(WeightFunctionArray=
zeros(length(ConvolvedDensity.DensityGridArray[1,:]), length(CosmologicalGrid.ZArray)),
BiasKind = CosmoCentral.PiecewiseBias())
CosmoCentral.ComputeBiasGrid!(CosmologicalGrid, GCWeightFunction, ConvolvedDensity)
x = CosmologicalGrid.ZArray; y = GCWeightFunction.BiasArray[1, :];
plot(x, y, ylabel = L"\mathrm{Bias}\,(z)", xlabel=L"z")
```
## Euclid Flagsip Bias
This model, obtained from the Flagship Euclid simulation is taken from
[Tutusaus et al 2021](https://arxiv.org/abs/2005.00055)
```math
b(z)=A+\frac{B}{1+\exp [-(z-D) C]}
```
```@example tutorial
GCWeightFunction = CosmoCentral.GCWeightFunction(WeightFunctionArray=
        zeros(length(ConvolvedDensity.DensityGridArray[1,:]),
        length(CosmologicalGrid.ZArray)), BiasKind = CosmoCentral.EuclidBias())
CosmoCentral.ComputeBiasGrid!(CosmologicalGrid, GCWeightFunction, ConvolvedDensity)
x = CosmologicalGrid.ZArray; y = GCWeightFunction.BiasArray[1, :];
plot(x, y, ylabel = L"\mathrm{Bias}\,(z)", xlabel= L"z")
```

```@docs
CosmoCentral.ComputeBias
CosmoCentral.ComputeBiasGrid!
```
