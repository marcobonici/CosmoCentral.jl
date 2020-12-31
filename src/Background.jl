"""
ComputeAdimensionalHubbleFactor(``z``, params)


This function, given the value of the cosmological parameters, evaluate the
Adimensional Hubble Factor for ``w_0 w_a``CDM cosmologies.

...
# Arguments
- `z::Float64` the redshift value at which evaluate the Adimensional Hubble Factor
- `params::w0waCDMCosmology`  a collection of cosmological parameters, whose expression is given by
```math
H(z)=H_0\\sqrt{\\Omega_M(1+z)^3+\\Omega_{DE}(1+z)^{3(1+w_0+w_a)\\exp(-3w_a \\frac{z}{1+z})}}
```


# Notes
* This expression is valid only for the CPL parameterization [LICENSE](https://github.com/marcobonici/CosmoCentral.jl/blob/main/LICENSE)
"""
function ComputeAdimensionalHubbleFactor(z::Float64, params::w0waCDMCosmology)
    E_z = sqrt(params.ΩM*(1+z)^3 + params.Ωr*(1+z)^4+params.Ωk*(1+z)^2
    +params.ΩDE*(1+z)^(3*(1+params.w0+params.wa))*exp(-3*params.wa*z/(1+z)) )
    return E_z
end

"""
ComputeHubbleFactor(``z``, params)

This function, given the value of the cosmological parameters, evaluate the
Hubble Factor for ``w_0 w_a``CDM cosmologies
"""
function ComputeHubbleFactor(z::Float64, params::w0waCDMCosmology)
    H_z = params.H0*ComputeAdimensionalHubbleFactor(z, params)
end

"""
ComputeComovingDistance(``z``, params)

This function, given the value of the cosmological parameters, evaluate the
Comoving Distance for ``w_0 w_a``CDM cosmologies. It is evaluated as
```math
r(z)=\\frac{c}{H_0}\\int_0^z \\frac{dx}{E(x)}
```
"""
function ComputeComovingDistance(z::Float64, params::w0waCDMCosmology)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    integral, err = QuadGK.quadgk(x -> 1 /
    ComputeAdimensionalHubbleFactor(x,params), 0, z, rtol=1e-12)
     return integral*c_0/params.H0
end
