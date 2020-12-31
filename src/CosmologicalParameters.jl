abstract type AbstractCosmology end
abstract type w0waCDMCosmology <: AbstractCosmology end

"""
w0waCDMParameters

This struct contains the value of the cosmological parameters for ``w_0 w_a``CDM cosmologies:
- ``w_0`` and ``w_a``, the parameters in the [CPL parameterization](https://arxiv.org/abs/astro-ph/0208512)

- ``\\Omega_M``, ``\\Omega_B``, ``\\Omega_{DE}``, ``\\Omega_R``, and ``\\Omega_k`` the density parameters for matter, baryons, Dark Energy, radiation, curvature

- ``n_s``, the scalar spectral index

- ``M_\\nu``, the sum of the neutrino mass eigenstates in eV

- ``\\sigma_8``, the amplitude of the scalar fluctuations

- ``H_0``, the value of the Hubble paramater
"""
functi
@kwdef struct w0waCDMParameters <: w0waCDMCosmology
    w0::Float64  = -1
    wa::Float64  = 0
    ΩM::Float64  = 0.32
    ΩB::Float64  = 0.05
    ΩDE::Float64 = 0.68
    Ωk::Float64  = 0.
    Ωr::Float64  = 0.
    ns::Float64  = 0.96
    Mν::Float64  = 0.06 #neutrino mass in eV
    σ8::Float64  = 0.816
    H0::Float64  = 67.
end
