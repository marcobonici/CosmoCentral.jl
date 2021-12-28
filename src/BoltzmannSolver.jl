"""
    Initializeclassy(cosmo::w0waCDMCosmology)

This function, given a [`w0waCDMCosmology`](@ref), returns the
[`classyParams`](@ref) correctly initialized.
"""
function Initializeclassy(cosmo::cosmo)
    classyParamsDict = Dict("output" => "mPk",
    "non linear"=> "halofit",
    "Omega_b"=> cosmo.ΩB,
    "Omega_cdm"=> cosmo.ΩM-cosmo.ΩB-
    cosmo.Mν/(93.14*(cosmo.H0/100)^2),
    "N_ur"=> 2.0328,
    "h"=> cosmo.H0/100.,
    "sigma8" => cosmo.σ8,
    "n_s" => cosmo.ns,
    "m_ncdm" => cosmo.Mν,
    "P_k_max_1/Mpc" => 50.,
    "z_max_pk" =>  5.,
    "use_ppf" =>  "yes",
    "w0_fld" =>  cosmo.w0,
    "Omega_k" =>  cosmo.Ωk,
    "Omega_fld" =>  cosmo.ΩDE,
    "wa_fld" =>  cosmo.wa,
    "cs2_fld" =>  1.,
    "N_ncdm" =>  1,
    "tau_reio" =>  0.058)
    classyparams = classyParams(classyParamsDict = classyParamsDict)
    return classyparams
end

function Initializeclassy(Flatw0waCDMCosmology::Flatw0waCDMCosmology)
    classyParamsDict = Dict("output" => "mPk",
    "non linear"=> "halofit",
    "Omega_b"=> Flatw0waCDMCosmology.ΩB,
    "Omega_cdm"=> Flatw0waCDMCosmology.ΩM-Flatw0waCDMCosmology.ΩB-
    Flatw0waCDMCosmology.Mν/(93.14*(Flatw0waCDMCosmology.H0/100)^2),
    "N_ur"=> 2.0328,
    "h"=> Flatw0waCDMCosmology.H0/100.,
    "sigma8" => Flatw0waCDMCosmology.σ8,
    "n_s" => Flatw0waCDMCosmology.ns,
    "m_ncdm" => Flatw0waCDMCosmology.Mν,
    "P_k_max_1/Mpc" => 50.,
    "z_max_pk" =>  5.,
    "use_ppf" =>  "yes",
    "w0_fld" =>  Flatw0waCDMCosmology.w0,
    "Omega_k" =>  0,
    "Omega_fld" =>  1-Flatw0waCDMCosmology.ΩM,
    "wa_fld" =>  Flatw0waCDMCosmology.wa,
    "cs2_fld" =>  1.,
    "N_ncdm" =>  1,
    "tau_reio" =>  0.058)
    classyparams = classyParams(classyParamsDict = classyParamsDict)
    return classyparams
end



"""
    EvaluatePowerSpectrum!(classyParams:: classyParams, cosmogrid::CosmologicalGrid,
    pmm::PowerSpectrum)

This function runs classy to evaluate the Matter Power Spectrum over the ``k-z``
grid specified in [`CosmologicalGrid`](@ref)
"""
function EvaluatePowerSpectrum!(classyParams:: classyParams,
    cosmogrid::CosmologicalGrid, pmm::PowerSpectrum)
    cosmo = classy.Class()
    cosmo.set(classyParams.classyParamsDict)
    cosmo.compute()
    for (idxz, myz) in enumerate(cosmogrid.ZArray)
        for (idxk, myk) in enumerate(cosmogrid.KArray)
            pmm.PowerSpectrumLinArray[idxk, idxz] =
            cosmo.pk_lin(myk, myz)
            pmm.PowerSpectrumNonlinArray[idxk, idxz] =
            cosmo.pk(myk, myz)
        end
    end
end
