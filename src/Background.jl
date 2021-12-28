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
function ComputeAdimensionalHubbleFactor(z::T, cosmo::w0waCDMCosmology) where T
    E_z = ComputeAdimensionalHubbleFactor(z, cosmo.ΩM, cosmo.Ωr,
    cosmo.ΩDE, cosmo.w0, cosmo.wa)
    return E_z
end

function ComputeAdimensionalHubbleFactor(z::T, ΩM::T, Ωr::T, ΩDE::T, w0::T, wa::T) where T
    Ωk = 1 - ΩM - Ωr - ΩDE
    E_z = sqrt(ΩM*(1+z)^3 + Ωr*(1+z)^4+ Ωk*(1+z)^2 
    +ΩDE*(1+z)^(3*(1+w0+wa))*exp(-3*wa*z/(1+z)) )
    return E_z
end

function ComputeAdimensionalHubbleFactor(z::T, 
    flatw0waCDMCosmology::Flatw0waCDMCosmology) where T
    ΩDE = 1 - cosmo.ΩM- cosmo.Ωr
    E_z = ComputeAdimensionalHubbleFactor(z, cosmo.ΩM, cosmo.Ωr,
    ΩDE, cosmo.w0, cosmo.wa)
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
function ComputeHubbleFactor(z::T, AbstractCosmology::AbstractCosmology) where T
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
function Computeχ(z::T, AbstractCosmology::AbstractCosmology) where T
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    integral, err = quadgk(x -> 1 /
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
        BackgroundQuantities.χZArray[idx_ZArray] = Computeχ(
        CosmologicalGrid.ZArray[idx_ZArray], AbstractCosmology)
    end
end

function ComputeLogSpacedχGrid!(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)
    Zχ = Dierckx.Spline1D(BackgroundQuantities.χZArray, CosmologicalGrid.ZArray)
    χArray = CosmoCentral.LogSpaced(BackgroundQuantities.χZArray[1],
    last(BackgroundQuantities.χZArray), length(BackgroundQuantities.χZArray))
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