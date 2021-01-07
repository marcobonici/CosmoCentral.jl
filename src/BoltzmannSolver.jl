

function EvaluatePowerSpectrum(classyParams:: classyParams)
    cosmo = classy.Class()
    cosmo.set(classyParams.classyParamsDict)
    cosmo.compute()
    for (idxz, myz) in enumerate(classyParams.CosmologicalGrid.ZArray)
        for (idxk, myk) in enumerate(classyParams.CosmologicalGrid.ZArray)
            classyParams.LinPowerSpectrumArray[idxk, idxz] =
            cosmo.pk_lin(myk, myz)
            classyParams.NonlinPowerSpectrumArray[idxk, idxz] =
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
