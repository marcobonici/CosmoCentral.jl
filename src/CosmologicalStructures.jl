abstract type AbstractCosmology{T} end
abstract type AbstractCosmologicalGrid{T} end
abstract type AbstractBackgroundQuantities{T} end
abstract type BoltzmannSolverParams end
abstract type AbstractDensity{T} end
abstract type AbstractConvolvedDensity{T} end
abstract type AbstractInstrumentResponse{T} end
abstract type AbstractWeightFunction{T} end
abstract type AbstractSourceFunction{T} end
abstract type AbstractTransferFunction end
abstract type AbstractPowerSpectrum{T} end
abstract type AbstractCℓ{T,N} end
abstract type Abstract∂Cℓ{T,N} end
abstract type AbstractBias{T} end
abstract type LensingEfficiencyMethod end
abstract type AbstractIntrinsicAlignment{T} end
abstract type AbstractFFTLog{T,C,I} end
abstract type AbstractCovariance{T} end
abstract type AbstractFisher{T,S} end
abstract type AbstractProbe end

"""
    w0waCDMCosmology(w0::Float64 = -1, wa::Float64 = 0, ΩM::Float64 = 0.32,
    ΩB::Float64  = 0.05, ΩDE::Float64 = 0.68, Ωk::Float64  = 0.,
    Ωr::Float64  = 0., ns::Float64  = 0.96, Mν::Float64  = 0.06,
    σ8::Float64  = 0.816, H0::Float64  = 67.)

This struct contains the value of the cosmological parameters for ``w_0 w_a``CDM cosmologies:
- ``w_0`` and ``w_a``, the parameters in the [CPL parameterization](https://arxiv.org/abs/astro-ph/0208512)

- ``\\Omega_M``, ``\\Omega_B``, ``\\Omega_{DE}``, ``\\Omega_R``, and ``\\Omega_k`` the density parameters for respectively  matter, baryons, Dark Energy, radiation, curvature

- ``n_s``, the scalar spectral index

- ``M_\\nu``, the sum of the neutrino mass eigenstates in eV

- ``\\sigma_8``, the amplitude of the scalar fluctuations

- ``H_0``, the value of the Hubble paramater
"""
@kwdef mutable struct w0waCDMCosmology{T} <: AbstractCosmology{T}
    #TODO #109
    w0::T  = -1.
    wa::T  = 0.
    Mν::T  = 0.06 #neutrino mass in eV
    H0::T  = 67.
    ΩM::T  = 0.32
    ΩB::T  = 0.05
    ΩDE::T = 0.68
    Ωk::T  = 0.
    Ωr::T  = 0.
    ns::T  = 0.96
    σ8::T  = 0.816
end

@kwdef mutable struct Flatw0waCDMCosmology{T} <: AbstractCosmology{T}
    w0::T  = -1.
    wa::T  = 0.
    Mν::T  = 0.06 #neutrino mass in eV
    H0::T  = 67.
    ΩM::T  = 0.32
    ΩB::T  = 0.05
    Ωr::T  = 0.
    ns::T  = 0.96
    σ8::T  = 0.816
end

"""
    CosmologicalGrid{T} <: AbstractCosmologicalGrid{T}
    ZArray::AbstractArray{T} = LinRange(0.001, 2.5, 300)
    KArray::AbstractArray{T} = LogSpaced(1e-5, 50., 1000)
    ℓBinCenters::AbstractArray{T} = LinRange(10., 3000., 2991)
    ℓBinWidths::AbstractArray{T} = LinRange(10., 3000., 2991)
    KLimberArray::AbstractArray{T,2} = zeros(length(ℓBinCenters),
    length(ZArray))
    KBeyondLimberArray::AbstractArray{T,2} = zeros(100, 1000)

This struct contains several grids:
- `ZArray` and `KArray`, the grids for the [`PowerSpectrum`](@ref) calculations
- `ℓBinCenters`, the ``\\ell`` grid for the ``C_\\ell``'s calculations
- `ℓBinWidths`, used when computing the covariance matrix
- `KLimberArray`, the `k` grid to evaluate the ``C_\\ell``'s in the Limber approximation
- `KBeyondLimberArray`, the `k` grid to evaluate the ``C_\\ell``'s without the Limber approximation
"""
@kwdef mutable struct CosmologicalGrid{T} <: AbstractCosmologicalGrid{T}
    ZArray::AbstractArray{T} = LinRange(0.001, 2.5, 300)
    KArray::AbstractArray{T} = LogSpaced(1e-5, 50., 1000)
    ℓBinCenters::AbstractArray{T} = LinRange(10., 3000., 2991)
    ℓBinWidths::AbstractArray{T} = LinRange(10., 3000., 2991)
    KLimberArray::AbstractArray{T,2} = zeros(length(ℓBinCenters),
    length(ZArray))
    KBeyondLimberArray::AbstractArray{T,2} = zeros(100, 1000)
end

"""
    BackgroundQuantities(HZArray::Vector{Float64}, χZArray::Vector{Float64}),
    DZArray::Vector{Float64}

This struct contains the arrays with the values of the Hubble parameter ``H(z)``
and the comoving distance ``\\chi(z)``.
"""
@kwdef mutable struct BackgroundQuantities{T} <: AbstractBackgroundQuantities{T}
    HZArray::Vector{T} = zeros(500)
    χZArray::Vector{T} = zeros(500)
    DZArray::Vector{T} = zeros(500)
end

"""
    classyParams(classyParamsDict::Dict)

This struct contains the dictionary with the classy parameters. For a detalied
explanation of the parameters, please refer to the
[CLASS website](http://class-code.net/)
"""
@kwdef mutable struct classyParams <: BoltzmannSolverParams
    classyParamsDict::Dict = Dict("output" => "mPk",
        "non linear"=> "halofit",
        "Omega_b"=> 0.05,
        "Omega_cdm"=> 0.2735,
        "N_ur"=> 2.0328,
        "h"=> 0.67,
        "sigma8" => 0.816,
        "n_s" => 0.96,
        "m_ncdm" => 0.06,
        "P_k_max_1/Mpc" => 50,
        "z_max_pk" =>  5,
        "use_ppf" =>  "yes",
        "w0_fld" =>  -1.,
        "Omega_k" =>  0,
        "Omega_fld" =>  0.68,
        "wa_fld" =>  0.,
        "cs2_fld" =>  1.,
        "N_ncdm" =>  1,
        "tau_reio" =>  0.058)
end

"""
    AnalitycalDensity(z0::Float64 = 0.9/sqrt(2.), zmin::Float64 = 0.001,
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

- `surfacedensity` , the value of the galaxy source density integrated between ``z_{min}`` and ``z_{max}``
- `normalization`, the value of parameter which multiplies the source dennsity in order to match the correct surface density
"""
@kwdef mutable struct AnalitycalDensity{T} <: AbstractDensity{T}
    Z0::T = 0.9/sqrt(2.)
    ZMin::T = 0.001
    ZMax::T = 4.0
    SurfaceDensity::T = 30.
    Normalization::T = 1.
end

"""
    ConvolvedDensity(AnalitycalDensity::AnalitycalDensity = AnalitycalDensity(),
    InstrumentResponse::InstrumentResponse = InstrumentResponse()
    ZBinArray::Vector{Float64} = Array([0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.50])
    DensityNormalizationArray::Vector{Float64} = ones(length(ZBinArray)-1)
    DensityGridArray::AbstractArray{Float64, 2} = ones(length(ZBinArray)-1, 300))


In order to take into account the error in the redshift measurement, the
source density is convolved with the [`InstrumentResponse`](@ref),
according to [the following equation](https://arxiv.org/abs/1910.09273)
```math
n_{i}(z)=\\frac{\\int_{z_{i}^{-}}^{z_{i}^{+}}
\\mathrm{d} z_{\\mathrm{p}} n(z) p \\left(z_{\\mathrm{p}}
\\mid z\\right)}{\\int_{z_{\\min }}^{z_{\\max }} \\mathrm{d} z
\\int_{z_{i}^{-}}^{z_{i}^{+}} \\mathrm{d} z_{\\mathrm{p}} n(z) p
\\left(z_{\\mathrm{p}} \\mid z\\right)}
```
"""
@kwdef mutable struct ConvolvedDensity{T} <: AbstractConvolvedDensity{T}
    ZBinArray::AbstractArray{T} = Array([0.001, 0.418, 0.560, 0.678, 0.789,
    0.900, 1.019, 1.155, 1.324, 1.576, 2.50])
    DensityNormalizationArray::AbstractArray{T} = ones(length(ZBinArray)-1)
    DensityGridArray::AbstractArray{T,2} = ones(length(ZBinArray)-1, 300)
    ShiftArray::AbstractArray{T} = zeros(length(ZBinArray)-1)
    SurfaceDensityArray::AbstractArray{T} = ones(10)
end

"""
    InstrumentResponse(cb::Float64 = 1.0, zb::Float64 = 0.0,
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
@kwdef struct InstrumentResponse{T} <: AbstractInstrumentResponse{T}
    cb::T   = 1.0
    zb::T   = 0.0
    σb::T   = 0.05
    co::T   = 1.0
    zo::T   = 0.1
    σo::T   = 0.05
    fout::T = 0.1
end

"""
    GCWeightFunction()

This struct contains the array with the Galaxy Bias and Galaxy Clustering Weight Function
values for all tomographic bins and redshift values in the [`CosmologicalGrid`](@ref).
"""
@kwdef mutable struct GCWeightFunction{T} <: AbstractWeightFunction{T}
    WeightFunctionArray::AbstractArray{T,2} = zeros(10, 500)
    BiasArray::AbstractArray{T,2} = zeros(size(WeightFunctionArray))
    BiasKind::AbstractBias = PiecewiseBias()
end


@kwdef mutable struct PiecewiseBias{T} <: AbstractBias{T}
    BiasMultiplier::AbstractArray{T} = ones(10)
end

"""
    EuclidBias()

This struct contains the parameter for the bias model measured by
[the Euclid Collaboration](https://arxiv.org/abs/2005.00055):

```math
b(z)= A +\\frac{B}{1+\\exp \\left( \\left(D-z \\right)C   \\right)}
```
"""
@kwdef mutable struct EuclidBias{T} <: AbstractBias{T}
    A::T = 1.0
    B::T = 2.5
    C::T = 2.8
    D::T = 1.6
end

struct AbsentIA{T} <: AbstractIntrinsicAlignment{T} end

@kwdef mutable struct ExtendedNLIA{T} <: AbstractIntrinsicAlignment{T}
    𝓐IA::T = 1.72
    βIA::T = 2.17
    𝓒IA::T = 0.0134
    ηIA::T = -0.41
end

"""
    WLWeightFunction()

This struct contains the array with the Lensing Efficiency and Weak Lensing
Weight Function values for all tomographic bins and redshift values in the
[`CosmologicalGrid`](@ref)
"""
@kwdef mutable struct WLWeightFunction{T} <: AbstractWeightFunction{T}
    WeightFunctionArray::AbstractArray{T,2} = zeros(10, 500)
    LensingEfficiencyArray::AbstractArray{T,2} = zeros(10, 500)
    IntrinsicAlignmentArray::AbstractArray{T,2} = zeros(10, 500)
    IntrinsicAlignmentModel::AbstractIntrinsicAlignment = ExtendedNLIA()
end

@kwdef mutable struct LensingSourceFunction{T} <: AbstractSourceFunction{T}
    SourceFunctionArray::AbstractArray{T,2} = zeros(10, 500)
    LensingEfficiencyArray::AbstractArray{T,2} = zeros(10, 500)
end


"""
    PowerSpectrum()

This struct contains the array with the Linear and Nonlinear Power Spectrum
evaluated on the ``k-z`` grid and the interpolated Nonlinear Power Spectrum on
Limber ``k-z`` grid.
"""
@kwdef mutable struct PowerSpectrum{T} <: AbstractPowerSpectrum{T}
    PowerSpectrumLinArray::AbstractArray{T,2} = zeros(1000, 300)
    PowerSpectrumNonlinArray::AbstractArray{T,2} = zeros(1000, 300)
    InterpolatedPowerSpectrum::AbstractArray{T,2} = zeros(2991, 300)
    GrowthFactor::AbstractArray{T} = zeros(
    length(PowerSpectrumLinArray[:,1]))
    InterpolatedPowerSpectrumBeyondLimber::AbstractArray{T,2} =
    zeros(100, 1000)
end

"""
    Cℓ(CℓArray::AbstractArray{Float64, 3})

This struct contains the array with the Angular Coefficients.
"""
@kwdef mutable struct Cℓ{T,N} <: AbstractCℓ{T,N}
    CℓArray::AbstractArray{T, N} = zeros(2991, 10, 10)
end

"""
    ∂Cℓ(∂CℓArray::AbstractArray{Float64, 3})

This struct contains the array with the derivatives of the Angular Coefficients.
"""
@kwdef mutable struct ∂Cℓ{T, N} <: Abstract∂Cℓ{T,N}
    ∂CℓArray::AbstractArray{T, N} = zeros(2991, 10, 10)
end

"""
    Fisherαβ()

This struct contains the array with the Fisher Matrix.
"""
@kwdef mutable struct Fisherαβ{T,S} <: AbstractFisher{T,S}
    FisherMatrix::AbstractArray{T,2} = zeros(8,8)
    FisherMatrixCumℓ::AbstractArray{T,3} = zeros(100,8,8)
    CorrelationMatrix::AbstractArray{T,2} = zeros(8,8)
    CorrelationMatrixCumℓ::AbstractArray{T,3} = zeros(100,8,8)
    FisherDict::Dict = Dict()
    FisherℓDict::Dict = Dict()
    CorrelationMatrixDict::Dict = Dict()
    ParametersList::Vector{S} = []
    SelectedParametersList::Vector{S} = []
    MarginalizedErrors::Dict = Dict()
    MarginalizedErrorsCumℓ::Dict = Dict()
    NormalizedErrors::Dict = Dict()
    NormalizedErrorsCumℓ::Dict = Dict()
end    

"""
    aₗₘCovariance()

This struct contains the Covariance in the field perspective, i.e. when the observables are the 
``a_{\\ell m}``'s.
"""
@kwdef mutable struct aₗₘCovariance{T}  <: AbstractCovariance{T}
    Covariance::AbstractArray{T,3} = zeros(2991, 10, 10)
    Cℓ::AbstractCℓ = Cℓ()
    Noise::AbstractArray{T,3} = zeros(2991, 10, 10)
    Covariance⁻¹::AbstractArray{T,3} = zeros(2991, 10, 10)
end

"""
    CℓCovariance()

This struct contains the Covariance in the estimator perspective, i.e. when the observables
are the ``C_{\\ell}``'s.
"""
@kwdef mutable struct CℓCovariance{T}  <: AbstractCovariance{T}
    Covariance::AbstractArray{T,3} = zeros(2991, 10, 10)
    Covariance⁻¹::AbstractArray{T,3} = zeros(2991, 10, 10)
end


"""
    InterpolationMethod

This type is used to specify the method to interpolate the Power Spectrum;
actually are included:

- Dierckx.jl, the Julia wrapper for the Fortran library Dierckx

- GriddedLinear, from Interpolations.jl

- BSpliceCubic (recommended, for its speed and accuracy), from Interpolations.jl
"""
abstract type InterpolationMethod end

struct RectBivSplineDierckx <: InterpolationMethod end
struct GriddedLinear <: InterpolationMethod end
struct BSplineCubic <: InterpolationMethod end

struct StandardLensingEfficiency <: LensingEfficiencyMethod end
struct CustomLensingEfficiency <: LensingEfficiencyMethod end

@kwdef mutable struct FFTLog{T,C,I} <: AbstractFFTLog{T,C,I}
    XArray::AbstractArray{T}
    DLnX::T = log(XArray[2]/XArray[1])
    FXArray::AbstractArray{T}
    OriginalLenght::I = length(XArray)
    ν::T = 1.01
    NExtrapLow::I = 0
    NExtrapHigh::I = 0
    CWindowWidth::T = 0.25
    NPad::I = 500
    N::I = OriginalLenght+NExtrapHigh+NExtrapLow+2*NPad
    M::AbstractArray{T} = zeros(N)
    CM::AbstractArray{C} = zeros(N)
    ηM::AbstractArray{T} = zeros(N)
end

"""
    κTransferFunction()

This struct contains the array with the Lensing Efficiency and Weak Lensing
Weight Function values for all tomographic bins and redshift values in the
[`CosmologicalGrid`](@ref)
"""
@kwdef mutable struct κTransferFunction <: AbstractTransferFunction
    LensingSourceFunction::LensingSourceFunction = LensingSourceFunction
    TransferFunctionArray::AbstractArray{Float64, 3} = zeros(10, 100, 1000)
end

@kwdef mutable struct GCProbe <: AbstractProbe
    Density::AbstractConvolvedDensity
    BiasModel::AbstractBias
    NuisanceDict::Dict = Dict()
end

@kwdef mutable struct WLProbe <:AbstractProbe
    Density::AbstractConvolvedDensity
    IAModel::AbstractIntrinsicAlignment
    NuisanceDict::Dict = Dict()
end