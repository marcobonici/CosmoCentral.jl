abstract type BackgroundQuantities end

"""
    ComputeAdimensionalHubbleFactor(z::Float64, params::w0waCDMCosmology)

This function, given the value of the cosmological parameters, evaluate the
Adimensional Hubble Factor for ``w_0 w_a``CDM cosmologies.
The analitycal expression is given by:
```math
E(z)=\\sqrt{\\Omega_M(1+z)^3+\\Omega_R(1+z)^4+\\Omega_{DE}(1+z)^{3(1+w_0+w_a)}\\exp\\left(-3w_a \\frac{z}{1+z}\\rigth)+\\Omega_k(1+z)^2}
```

!!! warning
    This expression is valid only for the
    [CPL parameterization](https://arxiv.org/abs/astro-ph/0208512)
    of the Dark Energy Equation of State.
"""
function ComputeAdimensionalHubbleFactor(z::Float64, cosmopars::w0waCDMCosmology)
    E_z = sqrt(cosmopars.ΩM*(1+z)^3 + cosmopars.Ωr*(1+z)^4+cosmopars.Ωk*(1+z)^2
    +cosmopars.ΩDE*(1+z)^(3*(1+cosmopars.w0+
    cosmopars.wa))*exp(-3*cosmopars.wa*z/(1+z)) )
    return E_z
end

"""
    ComputeHubbleFactor(z::Float64, params::w0waCDMCosmology)

This function, given the value of the cosmological parameters, evaluate the
Hubble Factor for ``w_0 w_a``CDM cosmologies, whose expression is given by
```math
H(z)=H_0\\sqrt{\\Omega_M(1+z)^3+\\Omega_R(1+z)^4+
\\Omega_{DE}(1+z)^{3(1+w_0+w_a)}\\exp(-3w_a \\frac{z}{1+z})+\\Omega_k(1+z)^2}
```

"""
function ComputeHubbleFactor(z::Float64, params::w0waCDMCosmology)
    H_z = params.H0*ComputeAdimensionalHubbleFactor(z, params)
end

"""
    ComputeComovingDistance(z::Float64, params::w0waCDMCosmology)

This function, given the value of the cosmological parameters, evaluate the
Comoving Distance. It is evaluated as:
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

"""
    BackgroundQuantitiesStruct(
    CosmologicalGrid::CosmologicalGrid = CosmologicalGridStruct()
    HZArray::Vector{Float64} = zeros(length(CosmologicalGrid.ZArray))
    rZArray::Vector{Float64} = zeros(length(CosmologicalGrid.ZArray))
    w0waCDMCosmology::w0waCDMCosmology = w0waCDMStruct())

This struct contains the value of the Cosmological Grid, both in ``k`` and ``z``.
"""
@kwdef struct BackgroundQuantitiesStruct <: BackgroundQuantities
    CosmologicalGrid::CosmologicalGrid = CosmologicalGridStruct()
    HZArray::Vector{Float64} = zeros(length(CosmologicalGrid.ZArray))
    rZArray::Vector{Float64} = zeros(length(CosmologicalGrid.ZArray))
    w0waCDMCosmology::w0waCDMCosmology = w0waCDMStruct()
end

function ComputeBackgroundQuantitiesOverGrid(
    BackgroundQuantities::BackgroundQuantities)
    for idx_ZArray in 1:length(BackgroundQuantities.CosmologicalGrid.ZArray)
        BackgroundQuantities.HZArray[idx_ZArray] = ComputeHubbleFactor(
        BackgroundQuantities.CosmologicalGrid.ZArray[idx_ZArray],
        BackgroundQuantities.w0waCDMCosmology)
        BackgroundQuantities.rZArray[idx_ZArray] = ComputeComovingDistance(
        BackgroundQuantities.CosmologicalGrid.ZArray[idx_ZArray],
        BackgroundQuantities.w0waCDMCosmology)
    end
end
