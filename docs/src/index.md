# CosmoCentral.jl

CosmoCentral is a Julia package to perform cosmological calculations. Actually
it can evaluate:

- Background quantities for ``w_0 w_a``CDM cosmologies
- Source densities, with an analitycal in input
- Angular Correlation functions, ``C_{ℓ}``'s, for several probes Galaxy
  Clustering

We aim to include also:
- Angular Correlation functions, ``C_{ℓ}``'s, for several probes (e.g., Weak
  Lensing, etc.)
- Fisher Matrix evaluation



### Authors

- Marco Bonici, Dipartimento di Fisica, Università degli Studi di Genova (UniGe)


## Usage

Here is an example of how CosmoCentral can be used to evaluate background
quantities.

```@repl
using CosmoCentral
w0waCDMCosmology = CosmoCentral.w0waCDMCosmologyStruct()
z = 1.
CosmoCentral.ComputeAdimensionalHubbleFactor(z, w0waCDMCosmology)
CosmoCentral.ComputeHubbleFactor(z, w0waCDMCosmology)
CosmoCentral.ComputeComovingDistance(z, w0waCDMCosmology)
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

### License

CosmoCentral is licensed under the MIT "Expat" license; see
[LICENSE](https://github.com/marcobonici/CosmoCentral.jl/blob/main/LICENSE) for
the full license text.
