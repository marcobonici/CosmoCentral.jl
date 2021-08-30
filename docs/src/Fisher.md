```@setup tutorial
using Plots; gr()
Plots.reset_defaults()
using CosmoCentral
```

# Fisher Forecast

Fisher Matrices can be used to forecast parameter uncertainties with few computational
resources. It is encapsulated in the following struct
```@docs
CosmoCentral.Fisherαβ
```
We implement two different approaches to evaluate the Fisher Matrix:

- A field approach, where we consider the ``a_{\ell m}`` to be the observable
- An estimator approach, where we consider the ``C_\ell``'s to be the observable

If correctly implemented, the two approaches gives the same result ([Hamimeche & Lewis 2008](https://arxiv.org/abs/0801.0554), [Carron 2012](https://arxiv.org/abs/1204.4724)).

## Field approach

The Fisher Matrix in the Field approach is evaluated by the following method

```@docs
CosmoCentral.ForecastFisherαβ(PathCentralCℓ::String, Path∂Cℓ::String,
InputList::Vector{Dict{String, Vector{Any}}},
CosmologicalGrid::CosmoCentral.CosmologicalGrid)
```

## Estimator approach

The Fisher Matrix in the Estimator approach is evaluated by the following method

```@docs
CosmoCentral.ForecastFisherαβ(PathCentralCℓ::String, Path∂Cℓ::String,
InputList::Vector{Dict{String, Vector{Any}}},
CosmologicalGrid::CosmoCentral.CosmologicalGrid, ciccio::String)
```

