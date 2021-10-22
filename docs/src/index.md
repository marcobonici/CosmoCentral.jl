# CosmoCentral.jl

CosmoCentral is a Julia package to perform cosmological calculations. Actually
it can evaluate:

- Background quantities for ``w_0 w_a``CDM cosmologies
- Galaxy redshift densities, with an analitycal in input
- Angular Correlation functions, ``C_{ℓ}``'s, for Galaxy Clustering and Weak Lensing using the Limber approximation
- Derivatives of ``C_{ℓ}``'s with respect to Cosmological Parameters
- Fisher Matrices and Marginalized errors

We aim to include also:
- Include other effects to the probes considered (e.g. magnification bias, redshift space distortions etc.)
- Evaluation of ``C_{ℓ}``'s beyond Limber approximation (WIP)
- Plot of Fisher Matrix 2D contours (WIP)



### Authors

- Marco Bonici, INAF - Institute of Space Astrophysics and Cosmic Physics (IASF), Milano
- Carmelita Carbone, INAF - Institute of Space Astrophysics and Cosmic Physics (IASF), Milano


## Usage

Here is an example of how CosmoCentral can be used to evaluate background
quantities.

```@repl
using CosmoCentral
w0waCDMCosmology = CosmoCentral.w0waCDMCosmology()
z = 1.
CosmoCentral.ComputeAdimensionalHubbleFactor(z, w0waCDMCosmology)
CosmoCentral.ComputeHubbleFactor(z, w0waCDMCosmology)
CosmoCentral.Computeχ(z, w0waCDMCosmology)
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

### License

CosmoCentral is licensed under the MIT "Expat" license; see
[LICENSE](https://github.com/marcobonici/CosmoCentral.jl/blob/main/LICENSE) for
the full license text.
