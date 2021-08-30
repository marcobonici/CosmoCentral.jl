"""
    ComputeAdimensionalHubbleFactor(z::Float64, w0waCDMCosmology::w0waCDMCosmology)

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
function ComputeAdimensionalHubbleFactor(z::Float64, w0waCDMCosmology::w0waCDMCosmology)
    E_z = sqrt(w0waCDMCosmology.ΩM*(1+z)^3 + w0waCDMCosmology.Ωr*(1+z)^4+
    w0waCDMCosmology.Ωk*(1+z)^2
    +w0waCDMCosmology.ΩDE*(1+z)^(3*(1+w0waCDMCosmology.w0+
    w0waCDMCosmology.wa))*exp(-3*w0waCDMCosmology.wa*z/(1+z)) )
    return E_z
end

function ComputeAdimensionalHubbleFactor(z::Float64,
    flatw0waCDMCosmology::Flatw0waCDMCosmology)
    E_z = sqrt(flatw0waCDMCosmology.ΩM*(1+z)^3
    +(1-flatw0waCDMCosmology.ΩM)*(1+z)^(3*(1+flatw0waCDMCosmology.w0+
    flatw0waCDMCosmology.wa))*exp(-3*flatw0waCDMCosmology.wa*z/(1+z)) )
    return E_z
end

"""
    ComputeHubbleFactor(z::Float64, AbstractCosmology::AbstractCosmology)

This function, given the value of the cosmological parameters, evaluate the
Hubble Factor for ``w_0 w_a``CDM cosmologies, whose expression is given by
```math
H(z)=H_0\\sqrt{\\Omega_M(1+z)^3+\\Omega_R(1+z)^4+
\\Omega_{DE}(1+z)^{3(1+w_0+w_a)}\\exp(-3w_a \\frac{z}{1+z})+\\Omega_k(1+z)^2}
```

"""
function ComputeHubbleFactor(z::Float64, AbstractCosmology::AbstractCosmology)
    H_z = AbstractCosmology.H0*ComputeAdimensionalHubbleFactor(z,
    AbstractCosmology)
    return H_z
end

"""
    Computeχ(z::Float64, AbstractCosmology::AbstractCosmology)

This function, given the value of the cosmological parameters, evaluate the
Comoving Distance χ. It is evaluated as:
```math
\\chi(z)=\\frac{c}{H_0}\\int_0^z \\frac{dz'}{E(z')}
```
"""
function Computeχ(z::Float64,
    AbstractCosmology::AbstractCosmology)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    integral, err = QuadGK.quadgk(x -> 1 /
    ComputeAdimensionalHubbleFactor(x,AbstractCosmology), 0, z, rtol=1e-12)
    return integral*c_0/AbstractCosmology.H0
end

"""
    ComputeBackgroundQuantitiesOverGrid(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    AbstractCosmology::AbstractCosmology)

This function evaluate the Hubble factor and the comoving distance over the
[`CosmologicalGrid`](@ref).
"""
function ComputeBackgroundQuantitiesGrid!(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    AbstractCosmology::AbstractCosmology)
    for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
        BackgroundQuantities.HZArray[idx_ZArray] = ComputeHubbleFactor(
        CosmologicalGrid.ZArray[idx_ZArray], AbstractCosmology)
        BackgroundQuantities.rZArray[idx_ZArray] = Computeχ(
        CosmologicalGrid.ZArray[idx_ZArray], AbstractCosmology)
    end
end

function ComputeLogSpacedχGrid!(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)
    Zχ = Dierckx.Spline1D(BackgroundQuantities.rZArray, CosmologicalGrid.ZArray)
    χArray = CosmoCentral.LogSpaced(BackgroundQuantities.rZArray[1],
    last(BackgroundQuantities.rZArray), length(BackgroundQuantities.rZArray))
    CosmologicalGrid.ZArray = Zχ.(χArray)
end

function ComputeBackgroundQuantitiesOverLogSpacedχGrid(
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    AbstractCosmology::AbstractCosmology)
    ComputeBackgroundQuantitiesGrid!(CosmologicalGrid, BackgroundQuantities,
    AbstractCosmology)
    ComputeLogSpacedχGrid!(CosmologicalGrid, BackgroundQuantities)
    ComputeBackgroundQuantitiesGrid!(CosmologicalGrid, BackgroundQuantities,
    AbstractCosmology)
end

function ExtractGrowthFactor!(BackgroundQuantities::BackgroundQuantities,
    PowerSpectrum::PowerSpectrum)
    BackgroundQuantities.DZArray = zeros(length(PowerSpectrum.PowerSpectrumLinArray[1,:]))
    for zidx in 1:length(BackgroundQuantities.DZArray)
        BackgroundQuantities.DZArray[zidx] =
        sqrt(PowerSpectrum.PowerSpectrumLinArray[1,zidx]./
        PowerSpectrum.PowerSpectrumLinArray[1,1])
    end
end