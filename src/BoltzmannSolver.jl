abstract type BoltzmannSolverParams end
abstract type classyParams <: BoltzmannSolverParams end

@kwdef mutable struct classyParamsStruct <: classyParams
    w0waCDMCosmology::w0waCDMCosmology = w0waCDMStruct()
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
        "tau_reio" =>  0.058,)
    PowerSpectrumGrid::PowerSpectrumGrid = PowerSpectrumGridStruct()
    LinPowerSpectrumArray::AbstractArray{Float64, 2} =
    zeros(length(PowerSpectrumGrid.kgrid), length(PowerSpectrumGrid.zgrid))
    NonlinPowerSpectrumArray::AbstractArray{Float64, 2} =
    zeros(length(PowerSpectrumGrid.kgrid), length(PowerSpectrumGrid.zgrid))
end

function EvaluatePowerSpectrum(classyParams:: classyParams)
    cosmo = classy.Class()
    cosmo.set(classyParams.classyParamsDict)
    cosmo.compute()
    for (idxz, myz) in enumerate(classyParams.PowerSpectrumGrid.zgrid)
        for (idxk, myk) in enumerate(classyParams.PowerSpectrumGrid.zgrid)
            classyParams.LinPowerSpectrumArray[idxk, idxz] =
            cosmo.pk_lin(myk, myz)
            classyParams.NonlinPowerSpectrumArray[idxk, idxz] =
            cosmo.pk(myk, myz)
        end
    end
end
