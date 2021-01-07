abstract type AbstractCosmology end
abstract type w0waCDMCosmology <: AbstractCosmology end
abstract type CosmologicalGrid end
abstract type AngularCoefficientsGrid end
abstract type BackgroundQuantities end
abstract type Bias end
abstract type PiecewiseBias <: Bias end
abstract type BoltzmannSolverParams end
abstract type classyParams <: BoltzmannSolverParams end
abstract type AsbtractDensity end
abstract type AnalitycalDensity <: AsbtractDensity end
abstract type ConvolvedDensity <: AsbtractDensity end
abstract type InstrumentResponse end
abstract type WeightFunction end
abstract type PowerSpectrum end




"""
    w0waCDMStruct(w0::Float64 = -1, wa::Float64 = 0, ΩM::Float64 = 0.32,
    ΩB::Float64  = 0.05, ΩDE::Float64 = 0.68, Ωk::Float64  = 0.,
    Ωr::Float64  = 0., ns::Float64  = 0.96, Mν::Float64  = 0.06,
    σ8::Float64  = 0.816, H0::Float64  = 67.)

This struct contains the value of the cosmological parameters for ``w_0 w_a``CDM cosmologies:
- ``w_0`` and ``w_a``, the parameters in the [CPL parameterization](https://arxiv.org/abs/astro-ph/0208512)

- ``\\Omega_M``, ``\\Omega_B``, ``\\Omega_{DE}``, ``\\Omega_R``, and ``\\Omega_k`` the density parameters for matter, baryons, Dark Energy, radiation, curvature

- ``n_s``, the scalar spectral index

- ``M_\\nu``, the sum of the neutrino mass eigenstates in eV

- ``\\sigma_8``, the amplitude of the scalar fluctuations

- ``H_0``, the value of the Hubble paramater
"""
@kwdef struct w0waCDMStruct <: w0waCDMCosmology
    w0::Float64  = -1.
    wa::Float64  = 0.
    Mν::Float64  = 0.06 #neutrino mass in eV
    H0::Float64  = 67.
    ΩM::Float64  = 0.32
    ΩB::Float64  = 0.05
    ΩDE::Float64 = 0.68
    Ωk::Float64  = 0.
    Ωr::Float64  = 0.
    ns::Float64  = 0.96
    σ8::Float64  = 0.816
end

"""
    CosmologicalGridStruct(ZArray::Vector{Float64} = Array(LinRange(0.001, 2.5, 300)),
    KArray::Vector{Float64} = LogSpaced(1e-5, 50., 1000))

This struct contains the value of the Cosmological Grid, both in ``k`` and ``z``.
"""
@kwdef struct CosmologicalGridStruct <: CosmologicalGrid
    ZArray::Vector{Float64} = Array(LinRange(0.001, 2.5, 300))
    KArray::Vector{Float64} = LogSpaced(1e-5, 50., 1000)
    MultipolesArray::Vector{Float64} = LinRange(10., 3000., 2991)
    KLimberArray::AbstractArray{Float64, 2} = zeros(length(MultipolesArray),
    length(ZArray))
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
end


@kwdef mutable struct PiecewiseBiasStruct <: PiecewiseBias
    BiasArray::AbstractArray{Float64, 2} =
    ones(10, 300)
end

@kwdef mutable struct classyParamsStruct <: classyParams
    classyParamsDict::Dict = Dict("output" => "mPk",
        "non linear"=> "halofit",
        "Omega_b"=> w0waCDMCosmology.ΩB,
        "Omega_cdm"=> w0waCDMCosmology.ΩM-w0waCDMCosmology.ΩB-
        w0waCDMCosmology.Mν/(93.14*(w0waCDMCosmology.H0/100)^2),
        "N_ur"=> 2.0328,
        "h"=> w0waCDMCosmology.H0/100.,
        "sigma8" => w0waCDMCosmology.σ8,
        "n_s" => w0waCDMCosmology.ns,
        "m_ncdm" => w0waCDMCosmology.Mν,
        "P_k_max_1/Mpc" => 50,
        "z_max_pk" =>  2.5,
        "use_ppf" =>  "yes",
        "w0_fld" =>  w0waCDMCosmology.w0,
        "Omega_k" =>  w0waCDMCosmology.Ωk,
        "Omega_fld" =>  w0waCDMCosmology.ΩDE,
        "wa_fld" =>  w0waCDMCosmology.wa,
        "cs2_fld" =>  1.,
        "N_ncdm" =>  1,
        "tau_reio" =>  0.058)
end

"""
    AnalitycalDensityStruct(z0::Float64 = 0.9/sqrt(2.), zmin::Float64 = 0.001,
    zmax::Float64 = 2.5, surfacedensity::Float64 = 30.,
    normalization::Float64 = 1.)

This struct contains the parameters of the source galaxy density as given by the
[official Euclid forecast](https://arxiv.org/abs/1910.09273), whose expression
is given by:
```math
n(z)\\propto\\left(\\frac{z}{z_0}\\right)^2
\\exp{\\left(-\\left(\\frac{z}{z_0}\\right)^{-3/2}\\right)}
```
The parameters contained in this struct are
- ``z_{min}`` and ``z_{max}``, the minimum and the maximum redshift considered

- ``z_0``, the parameter present in the galaxy distribution

- surfacedensity , the value of the galaxy source density integrated between ``z_{min}`` and ``z_{max}``
- normalization, the value of parameter which multiplies the source dennsity in order to match the correct surface density
"""
@kwdef mutable struct AnalitycalDensityStruct <: AnalitycalDensity
    Z0::Float64 = 0.9/sqrt(2.)
    ZMin::Float64 = 0.001
    ZMax::Float64 = 2.5
    SurfaceDensity::Float64 = 30.
    Normalization::Float64 = 1.
end

"""
    ConvolvedDensityStruct(AnalitycalDensity::AnalitycalDensity = AnalitycalDensityStruct(),
    InstrumentResponse::InstrumentResponse = InstrumentResponseStruct()
    ZBinArray::Vector{Float64} = Array([0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.50])
    DensityNormalizationArray::Vector{Float64} = ones(length(ZBinArray)-1)
    DensityGridArray::AbstractArray{Float64, 2} = ones(length(ZBinArray)-1, 300))


In order to take into account the error in the redshift measurement, the
source density is convolved with the [`InstrumentResponseStruct`](@ref),
according to [the following equation](https://arxiv.org/abs/1910.09273)
```math
n_{i}(z)=\\frac{\\int_{z_{i}^{-}}^{z_{i}^{+}}
\\mathrm{d} z_{\\mathrm{p}} n(z) p \\left(z_{\\mathrm{p}}
\\mid z\\right)}{\\int_{z_{\\min }}^{z_{\\max }} \\mathrm{d} z
\\int_{z_{i}^{-}}^{z_{i}^{+}} \\mathrm{d} z_{\\mathrm{p}} n(z) p
\\left(z_{\\mathrm{p}} \\mid z\\right)}
```
"""
@kwdef mutable struct ConvolvedDensityStruct <: ConvolvedDensity
    ZBinArray::Vector{Float64} = Array([0.001, 0.418, 0.560, 0.678, 0.789,
    0.900, 1.019, 1.155, 1.324, 1.576, 2.50])
    DensityNormalizationArray::Vector{Float64} = ones(length(ZBinArray)-1)
    DensityGridArray::AbstractArray{Float64, 2} = ones(length(ZBinArray)-1, 300)
end

"""
    InstrumentResponseStruct(cb::Float64 = 1.0, zb::Float64 = 0.0,
    σb::Float64 = 0.05, co::Float64 = 1.0, zo::Float64 = 0.1,
    σo::Float64 = 0.05, fout::Float64 = 0.1)

When we measure the redshift of a galaxy with redshit ``z``, we will measure a
redshift ``z_p`` with a probability given by
[the following expression](https://arxiv.org/abs/1910.09273):
```math
p(z_p|z)  = \\frac{1-f_{out}}{\\sqrt{2 \\pi} \\sigma_{b}(1+z)} \\exp \\left(
-\\frac{1}{2}\\left(\\frac{z-c_{b} z_{b}-z_{b}}{\\sigma_{b}(1+z)}\\right)^{2}
\\right) + \\frac{f_{out}}{\\sqrt{2 \\pi} \\sigma_{\\mathrm{o}}(1+z)} \\exp
\\left(-\\frac{1}{2}\\left(\\frac{z-c_{o} z_{p}-z_{o}}{\\sigma_{o}(1+z)}
\\right)^{2}\\right)
```
This struct contains all these parameters.
"""
@kwdef struct InstrumentResponseStruct <: InstrumentResponse
    cb::Float64   = 1.0
    zb::Float64   = 0.0
    σb::Float64   = 0.05
    co::Float64   = 1.0
    zo::Float64   = 0.1
    σo::Float64   = 0.05
    fout::Float64 = 0.1
end

@kwdef mutable struct WeightFunctionStruct <: WeightFunction
    ConvolvedDensity::ConvolvedDensity = ConvolvedDensityStruct()
    CosmologicalGrid::CosmologicalGrid = CosmologicalGridStruct()
    w0waCDMCosmology::w0waCDMCosmology = w0waCDMStruct()
    Bias::Bias = PiecewiseBiasStruct()
    BackgroundQuantities::BackgroundQuantities = BackgroundQuantitiesStruct()
    WeightFunctionArray::AbstractArray{Float64, 2} =
    ones(length(ConvolvedDensity.ZBinArray)-1, length(CosmologicalGrid.ZArray))
end

@kwdef mutable struct PowerSpectrumStruct <: PowerSpectrum
    PowerSpectrumLinArray::AbstractArray{Float64, 2} =
    zeros(1000, 300)
    PowerSpectrumNonlinArray::AbstractArray{Float64, 2} =
    zeros(1000, 300)
end