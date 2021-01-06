abstract type AbstractCosmology end
abstract type w0waCDMCosmology <: AbstractCosmology end
abstract type CosmologicalGrid end
abstract type AngularCoefficientsGrid end
abstract type BackgroundQuantities end



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

function ComputeKLimberArray(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)
    for idx_z in 1:length((CosmologicalGrid.ZArray))
        CosmologicalGrid.KLimberArray[:, myz] =
        (CosmologicalGrid.MultipolesArray.+
        1. /2.)./BackgroundQuantities.HZArray[idx_z]
    end
end
