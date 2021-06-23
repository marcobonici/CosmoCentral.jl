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

In real surveys we do not deal with the exactr distributions due to errors in
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
