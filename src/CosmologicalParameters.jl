abstract type AbstractCosmology end
abstract type w0waCDMCosmology <: AbstractCosmology end


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
