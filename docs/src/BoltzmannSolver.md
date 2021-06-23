# BoltzmannSolver

One of the core components of CosmoCentral are the Boltzmann solvers, codes
developed by the community that solve the coupled Einstein-Boltzmann equations
and used to evaluate the Matter Power Spectrum. Currently, the Boltzmann codes
implemented in CosmoCentral are:

- CLASS

We plan to include other Boltzmann solvers in the future.

## Power Spectrum

```@docs
CosmoCentral.PowerSpectrum
CosmoCentral.ComputeLimberArray!
CosmoCentral.InterpolatePowerSpectrumLimberGrid!
CosmoCentral.InterpolationMethod
```

## classy

classy is the Python wrapper for [CLASS](http://class-code.net/).

```@docs
CosmoCentral.classyParams
CosmoCentral.Initializeclassy
CosmoCentral.EvaluatePowerSpectrum!
```
