"""
ComputeAdimensionalHubbleFactor(z, ``w_0``, ``wa``, ``\\Omega_M``, ``\\Omega_{DE}``, ``\\Omega_k``, ``\\Omega_R``)


This function, given the value of the cosmological parameters, evaluate the
Adimensional Hubble Factor for w0waCDM cosmologies

...
# Arguments
- `n::Integer`: the number of elements to compute.
- `dim::Integer=1`: the dimensions along which to perform the computation.
...

# Notes
* Notes can go here

# Examples
```julia
julia> five = plusTwo(3)
5

"""
function ComputeAdimensionalHubbleFactor(z::Float64, w0::Float64=-1.0,
    wa::Float64=0.0, ΩM::Float64=0.32, ΩDE::Float64=0.68, Ωk::Float64=0.0,
    Ωr::Float64=0.0)
    E_z = sqrt(ΩM*(1+z)^3 + Ωr*(1+z)^4+Ωk*(1+z)^2
    +ΩDE*(1+z)^(3*(1+w0+wa))*exp(-3*wa*z/(1+z)) )
    return E_z
end

"""
ComputeHubbleFactor(z, H0, w0, wa, ΩM, ΩDE, Ωk, Ωr)

This function, given the value of the cosmological parameters, evaluate the
Hubble Factor for w0waCDM cosmologies
"""
function ComputeHubbleFactor(z::Float64, H0::Float64=67.0, w0::Float64=-1.0,
    wa::Float64=0.0, ΩM::Float64=0.32, ΩDE::Float64=0.68, Ωk::Float64=0.0,
    Ωr::Float64=0.0)
    H_z = H0*ComputeAdimensionalHubbleFactor(z, w0, wa, ΩM, ΩDE, Ωk, Ωr)
end

"""
ComputeComovingDistance(z, H0, w0, wa, ΩM, ΩDE, Ωk, Ωr)

This function, given the value of the cosmological parameters, evaluate the
Comoving Distance for w0waCDM cosmologies
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
