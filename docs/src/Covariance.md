```@setup tutorial
using Plots; gr()
Plots.reset_defaults()
using CosmoCentral
```

# Covariance Matrix

In order to evaluate the Fisher Matrix [`CosmoCentral.Fisherαβ`](@ref), a Covariance Matrix
is required. Here we show the two Gaussian Covariance Matrix we have implemented, in the
Field and Estimator Approach.

## Field approach

The Covariance Matrix in the Field approach is evaluated by the following method

```@docs
CosmoCentral.aₗₘCovariance
CosmoCentral.InstantiateEvaluateCovariance(cℓ::CosmoCentral.AbstractCℓ,
ConvDens::CosmoCentral.AbstractConvolvedDensity, cosmogrid::CosmoCentral.CosmologicalGrid,
ProbeA::String, ProbeB::String)
```

## Estimator approach

The Fisher Matrix in the Estimator approach is evaluated by the following method

```@docs
CosmoCentral.CℓCovariance
CosmoCentral.InstantiateEvaluateCovariance(Covaₗₘ::CosmoCentral.aₗₘCovariance)
```