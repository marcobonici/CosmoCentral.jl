# Cosmology
In this section we show the structures which specify the cosmology and the value of the
parameters. In this moment we support only the ``w_0 w_a``CDM cosmology, although we plan to
include support to other cosmologies.

```@docs
CosmoCentral.w0waCDMCosmology
```
In the remainder of this documentation, we will show how to use CosmoCentral.jl and we will
plot several quantities. The assumed values for the cosmological parameters are:
- ``\Omega_M`` = 0.32
- ``\Omega_B`` = 0.05
- ``\Omega_k`` = 0
- ``n_s`` = 0.96
- ``\sigma_8`` = 0.816
- ``h`` = 0.67
- ``M_\nu`` = 0.06 eV
- ``w_0`` = -1
- ``w_a`` = 0

## Cosmological Grid
Another import quantity is the cosmological grid, the structure which specifies the grids
used to compute the quantities in our calculations.
```@docs
CosmoCentral.CosmologicalGrid
```