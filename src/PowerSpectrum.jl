"""
    ComputeLimberArray(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)

This function compute the Limber grid. In the Limber approximation, ``k`` and
``z`` are related by the following relation:
```math
k_\\ell(z)=\\frac{\\ell+1/2}{r(z)}.
```
"""
function ComputeLimberArray(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)
    for (idx_l, myl) in enumerate(CosmologicalGrid.MultipolesArray)
        CosmologicalGrid.KLimberArray[idx_l, :] =  (myl + 0.5) ./
        BackgroundQuantities.rZArray
    end
end

"""
    InterpolateAndEvaluatePowerSpectrum(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, PowerSpectrum::PowerSpectrum)

This function interpolates the Power Spectrum on the ``k-z`` grid and evaluates
it on the Limber grid.
"""
function InterpolateAndEvaluatePowerSpectrum(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, PowerSpectrum::PowerSpectrum,
    ::RectBivSplineDierckx)
    InterpPmm = Dierckx.Spline2D(log10.(CosmologicalGrid.KArray),
    CosmologicalGrid.ZArray, log10.(PowerSpectrum.PowerSpectrumNonlinArray);
    kx=5, ky=5, s=0.0)
    for (idx_l, myl) in enumerate(CosmologicalGrid.MultipolesArray)
        PowerSpectrum.InterpolatedPowerSpectrum[idx_l, :] =
        10 .^(InterpPmm.(log10.(CosmologicalGrid.KLimberArray[idx_l, :]),
        CosmologicalGrid.ZArray))
    end
end

function InterpolateAndEvaluatePowerSpectrum(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, PowerSpectrum::PowerSpectrum,
    ::GriddedLinear)
    InterpPmm = Interpolations.interpolate((log10.(CosmologicalGrid.KArray),
    CosmologicalGrid.ZArray), log10.(PowerSpectrum.PowerSpectrumNonlinArray),
    Gridded(Linear()))
    InterpPmm = Interpolations.extrapolate(InterpPmm, Line())
    for (idx_l, myl) in enumerate(CosmologicalGrid.MultipolesArray)
        PowerSpectrum.InterpolatedPowerSpectrum[idx_l, :] =
        10 .^(InterpPmm.(log10.(CosmologicalGrid.KLimberArray[idx_l, :]),
        CosmologicalGrid.ZArray))
    end
end
