<p align="center">
<img width="400px" src="https://raw.githubusercontent.com/marcobonici/CosmoCentral.jl/develop/docs/src/assets/logo.png"/>
</p>

# CosmoCentral.jl

CosmoCentral is a code for cosmological analysis written in Julia.
| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://marcobonici.github.io/CosmoCentral.jl/dev) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://marcobonici.github.io/CosmoCentral.jl/stable) | [![Build status (Github Actions)](https://github.com/marcobonici/CosmoCentral.jl/workflows/CI/badge.svg)](https://github.com/marcobonici/CosmoCentral.jl/actions) [![codecov](https://codecov.io/gh/marcobonici/CosmoCentral.jl/branch/main/graph/badge.svg?token=VK3FKFHMVQ)](https://codecov.io/gh/marcobonici/CosmoCentral.jl)|

## Usage

```python
import CosmoCentral

params = CosmoCentral.w0waCDMParameters()

CosmoCentral.ComputeAdimensionalHubbleFactor(z, params) # returns the adimensional Hubble factor
CosmoCentral.ComputeHubbleFactor(z, params) # returns the Hubble factor
CosmoCentral.ComputeComovingDistance(z, params) # returns the comoving distance
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.
