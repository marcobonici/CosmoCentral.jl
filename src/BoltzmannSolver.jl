function Initializeclassy(w0waCDMCosmology::w0waCDMCosmology)
    classyParamsDict = Dict("output" => "mPk",
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
    classyParams = classyParamsStruct(classyParamsDict = classyParamsDict)
end


function EvaluatePowerSpectrum(classyParams:: classyParams,
    CosmologicalGrid::CosmologicalGrid, PowerSpectrum::PowerSpectrum)
    cosmo = classy.Class()
    cosmo.set(classyParams.classyParamsDict)
    cosmo.compute()
    for (idxz, myz) in enumerate(CosmologicalGrid.ZArray)
        for (idxk, myk) in enumerate(CosmologicalGrid.KArray)
            PowerSpectrum.PowerSpectrumLinArray[idxk, idxz] =
            cosmo.pk_lin(myk, myz)
            PowerSpectrum.PowerSpectrumNonlinArray[idxk, idxz] =
            cosmo.pk(myk, myz)
        end
    end
end

function ComputeKLimberArray(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)
    for idx_z in 1:length((CosmologicalGrid.ZArray))
        CosmologicalGrid.KLimberArray[:, myz] =
        (CosmologicalGrid.MultipolesArray.+
        1. /2.)./BackgroundQuantities.HZArray[idx_z]
    end
end
