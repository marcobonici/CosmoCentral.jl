"""
ComputeAdimensionalHubbleFactor(``z``, ``w_0``, ``w_a``, ``\\Omega_M``, ``\\Omega_{DE}``, ``\\Omega_k``, ``\\Omega_R``)


This function, given the value of the cosmological parameters, evaluate the
Adimensional Hubble Factor for ``w_0 w_a``CDM cosmologies

...
# Arguments
- `z::Float64`: the redshift value at which evaluate the Adimensional Hubble Factor


# Notes
* This expression is valid only for the CPL parameterization

# Examples
```julia
julia> ComputeAdimensionalHubbleFactor(1.0)
1.0

"""
function ComputeAdimensionalHubbleFactor(z::Float64, w0::Float64=-1.0,
    wa::Float64=0.0, ΩM::Float64=0.32, ΩDE::Float64=0.68, Ωk::Float64=0.0,
    Ωr::Float64=0.0)
    E_z = sqrt(ΩM*(1+z)^3 + Ωr*(1+z)^4+Ωk*(1+z)^2
    +ΩDE*(1+z)^(3*(1+w0+wa))*exp(-3*wa*z/(1+z)) )
    return E_z
end

"""
ComputeHubbleFactor(``z``, ``H_0``,``w_0``, ``w_a``, ``\\Omega_M``, ``\\Omega_{DE}``, ``\\Omega_k``, ``\\Omega_R``)

This function, given the value of the cosmological parameters, evaluate the
Hubble Factor for ``w_0 w_a``CDM cosmologies
"""
function ComputeHubbleFactor(z::Float64, H0::Float64=67.0, w0::Float64=-1.0,
    wa::Float64=0.0, ΩM::Float64=0.32, ΩDE::Float64=0.68, Ωk::Float64=0.0,
    Ωr::Float64=0.0)
    H_z = H0*ComputeAdimensionalHubbleFactor(z, w0, wa, ΩM, ΩDE, Ωk, Ωr)
end

"""
ComputeComovingDistance(``z``, ``H_0``,``w_0``, ``w_a``, ``\\Omega_M``, ``\\Omega_{DE}``, ``\\Omega_k``, ``\\Omega_R``)

This function, given the value of the cosmological parameters, evaluate the
Comoving Distance for ``w_0 w_a``CDM cosmologies. It is evaluated as
```math
r(z)=\\frac{c}{H_0}\\int_0^z \\frac{dx}{E(x)}
```
"""
function ComputeComovingDistance(z::Float64, H0::Float64=67.0, w0::Float64=-1.0,
    wa::Float64=0.0, ΩM::Float64=0.32, ΩDE::Float64=0.68, Ωk::Float64=0.0,
    Ωr::Float64=0.0,)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    integral, err = QuadGK.quadgk(x -> 1 / ComputeAdimensionalHubbleFactor(x,
    w0, wa, ΩM, ΩDE, Ωk, Ωr), 0, z, rtol=1e-12)

     return integral*c_0/H0
end
