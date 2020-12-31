# CosmoCentral.jl

Documentation for CosmoCentral.jl. This is a Julia package to perform cosmological calculations.



### Authors

- Marco Bonici, Dipartimento di Fisica, Università degli Studi di Genova (UniGe)


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

### License

CosmoCentral is licensed under the MIT "Expat" license; see
[LICENSE](https://github.com/marcobonici/CosmoCentral.jl/blob/main/LICENSE) for
the full license text.
