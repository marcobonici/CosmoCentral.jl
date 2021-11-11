```@setup tutorial
using Plots; gr()
Plots.reset_defaults()
using PlotThemes
using CosmoCentral
using LaTeXStrings

#instantiate cosmology and grid
w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()
CosmologicalGrid  = CosmoCentral.CosmologicalGrid(
ZArray=Array(LinRange(0.001, 4.0, 500)))

plot_font = "Computer Modern"
Plots.default(titlefont = (16, plot_font), fontfamily=plot_font,
        linewidth=2, framestyle=:box, fg_legend =:black, label=nothing, grid=false,
        tickfontsize=12, legendfontsize=12, size = (550, 400), labelfontsize = 13,
        dpi = 200)

```


# Source Density
In this page are presented the structures and functions used to deal with source
densities and are listed here.

```@index
Pages = ["SourceDensity.md"]
```

## Source Density

In the evaluation of Angular Coefficients, central quantities are the source
densities. In this section are presented the custom types and function used to
deal with the source densities.

```@docs
CosmoCentral.AnalitycalDensity
CosmoCentral.ComputeDensity
CosmoCentral.NormalizeAnalitycalDensity!
```

## Convolved Source Density

In real surveys we do not deal with the exact distributions due to errors in
the measurement of the source redshifts. The redshift errors are accounted for
convolving the source density with a redshift measurement error.

### Intrument Response

```@docs
CosmoCentral.InstrumentResponse
CosmoCentral.ComputeInstrumentResponse
```

### Convolved Source Density

```@docs
CosmoCentral.ConvolvedDensity
CosmoCentral.ComputeConvolvedDensity
CosmoCentral.NormalizeConvolvedDensity!
CosmoCentral.ComputeConvolvedDensityGrid!
```
Here we show how to calculate ``n_g^i(z)``, then we will plot it
```@example tutorial

#instantiate the analytical density and normalize it
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

p = Plots.plot(xlabel=L"z", ylabel=L"n_i^g(z)",
    title="Normalized galaxy density")
for i in 1:10
Plots.plot!(p, CosmologicalGrid.ZArray, ConvolvedDensity.DensityGridArray[i,:],
    labels=(L"i=%$i"),  linewidth=3)
end
p
```