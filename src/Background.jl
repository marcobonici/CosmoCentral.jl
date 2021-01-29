"""
    ComputeAdimensionalHubbleFactor(z::Float64, params::w0waCDMCosmologyStruct)

This function, given the value of the cosmological parameters, evaluate the
Adimensional Hubble Factor for ``w_0 w_a``CDM cosmologies.
The analitycal expression is given by:
```math
E(z)=\\sqrt{\\Omega_M(1+z)^3+\\Omega_R(1+z)^4+
\\Omega_{DE}(1+z)^{3(1+w_0+w_a)}\\exp(-3w_a \\frac{z}{1+z})+\\Omega_k(1+z)^2}
```

!!! warning
    This expression is valid only for the
    [CPL parameterization](https://arxiv.org/abs/astro-ph/0208512)
    of the Dark Energy Equation of State.
"""
function ComputeAdimensionalHubbleFactor(z::Float64,
    w0waCDMCosmology::w0waCDMCosmologyStruct)
    E_z = sqrt(w0waCDMCosmology.ΩM*(1+z)^3 + w0waCDMCosmology.Ωr*(1+z)^4+
    w0waCDMCosmology.Ωk*(1+z)^2
    +w0waCDMCosmology.ΩDE*(1+z)^(3*(1+w0waCDMCosmology.w0+
    w0waCDMCosmology.wa))*exp(-3*w0waCDMCosmology.wa*z/(1+z)) )
    return E_z
end

"""
    ComputeHubbleFactor(z::Float64, params::w0waCDMCosmologyStruct)

This function, given the value of the cosmological parameters, evaluate the
Hubble Factor for ``w_0 w_a``CDM cosmologies, whose expression is given by
```math
H(z)=H_0\\sqrt{\\Omega_M(1+z)^3+\\Omega_R(1+z)^4+
\\Omega_{DE}(1+z)^{3(1+w_0+w_a)}\\exp(-3w_a \\frac{z}{1+z})+\\Omega_k(1+z)^2}
```

"""
function ComputeHubbleFactor(z::Float64,
    w0waCDMCosmology::w0waCDMCosmologyStruct)
    H_z = w0waCDMCosmology.H0*ComputeAdimensionalHubbleFactor(z,
    w0waCDMCosmology)
    return H_z
end

"""
    ComputeComovingDistance(z::Float64, params::w0waCDMCosmology)

This function, given the value of the cosmological parameters, evaluate the
Comoving Distance. It is evaluated as:
```math
r(z)=\\frac{c}{H_0}\\int_0^z \\frac{dz'}{E(z')}
```
"""
function ComputeComovingDistance(z::Float64,
    w0waCDMCosmology::w0waCDMCosmologyStruct)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    integral, err = QuadGK.quadgk(x -> 1 /
    ComputeAdimensionalHubbleFactor(x,w0waCDMCosmology), 0, z, rtol=1e-12)
    return integral*c_0/w0waCDMCosmology.H0
end

"""
    ComputeBackgroundQuantitiesOverGrid(CosmologicalGrid::CosmologicalGrid,
        BackgroundQuantities::BackgroundQuantities,
        w0waCDMCosmology::w0waCDMCosmology)

This function evaluate the Hubble factor and the comoving distance over the
[`CosmologicalGridStruct`](@ref).
"""
function ComputeBackgroundQuantitiesOverGrid(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct)
    for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
        BackgroundQuantities.HZArray[idx_ZArray] = ComputeHubbleFactor(
        CosmologicalGrid.ZArray[idx_ZArray], w0waCDMCosmology)
        BackgroundQuantities.rZArray[idx_ZArray] = ComputeComovingDistance(
        CosmologicalGrid.ZArray[idx_ZArray], w0waCDMCosmology)
    end
end
