# BoltzmannSolver

One of the core components of CosmoCentral are the Boltzmann solvers, codes
developed by the community that solve the coupled Einstein-Boltzmann equations
and used to evaluate the Matter Power Spectrum. Currently, the Boltzmann codes
implemented in CosmoCentral are:

- CLASS

We plan to include other Boltzmann solvers in the future.

## Power Spectrum

```@docs
CosmoCentral.PowerSpectrumStruct
CosmoCentral.ComputeLimberArray
CosmoCentral.InterpolateAndEvaluatePowerSpectrum
CosmoCentral.InterpolationMethod
```

## classy

classy is the Python wrapper for [CLASS](http://class-code.net/).

```@docs
CosmoCentral.classyParamsStruct
CosmoCentral.Initializeclassy
CosmoCentral.EvaluatePowerSpectrum
```