"""
    ComputeLimberArray!(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)

This function compute the Limber grid. In the Limber approximation, ``k`` and
``z`` are related by the following relation:
```math
k_\\ell(z)=\\frac{\\ell+1/2}{r(z)}.
```
"""
function ComputeLimberArray!(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)
    for (idx_l, myl) in enumerate(CosmologicalGrid.ℓBinCenters)
        CosmologicalGrid.KLimberArray[idx_l, :] =  (myl + 0.5) ./
        BackgroundQuantities.χZArray
    end
end

"""
    InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, PowerSpectrum::PowerSpectrum,
    ::InterpolationMethod)

This function interpolates the Power Spectrum on the ``k-z`` grid and evaluates
it on the Limber grid.
"""
function InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, PowerSpectrum::PowerSpectrum,
    ::RectBivSplineDierckx)
    InterpPmm = Dierckx.Spline2D(log10.(CosmologicalGrid.KArray),
    CosmologicalGrid.ZArray, log10.(PowerSpectrum.PowerSpectrumNonlinArray);
    kx=5, ky=5, s=0.0)
    @inbounds @simd for idx_l in 1:length(CosmologicalGrid.ℓBinCenters)
        PowerSpectrum.InterpolatedPowerSpectrum[idx_l, :] =
        10 .^(InterpPmm.(log10.(CosmologicalGrid.KLimberArray[idx_l, :]),
        CosmologicalGrid.ZArray))
    end
end

function InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, PowerSpectrum::PowerSpectrum,
    ::BSplineCubic)
    x = LinRange(log10(first(CosmologicalGrid.KArray)),
    log10(last(CosmologicalGrid.KArray)), length(CosmologicalGrid.KArray))
    y = LinRange(first(CosmologicalGrid.ZArray), last(CosmologicalGrid.ZArray),
    length(CosmologicalGrid.ZArray))
    InterpPmm = Interpolations.interpolate(
    log10.(PowerSpectrum.PowerSpectrumNonlinArray),
    BSpline(Cubic(Line(OnGrid()))))
    InterpPmm = scale(InterpPmm, x, y)
    InterpPmm = Interpolations.extrapolate(InterpPmm, Line())
    @inbounds @simd for idx_l in 1:length(CosmologicalGrid.ℓBinCenters)
        PowerSpectrum.InterpolatedPowerSpectrum[idx_l, :] =
        10 .^(InterpPmm.(log10.(CosmologicalGrid.KLimberArray[idx_l, :]),
        CosmologicalGrid.ZArray))
    end
end

function InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, PowerSpectrum::PowerSpectrum,
    ConvolvedDensity::ConvolvedDensity)
    x = LinRange(log10(first(CosmologicalGrid.KArray)),
    log10(last(CosmologicalGrid.KArray)), length(CosmologicalGrid.KArray))
    y = LinRange(first(CosmologicalGrid.ZArray), last(CosmologicalGrid.ZArray),
    length(CosmologicalGrid.ZArray))
    InterpPmm = Interpolations.interpolate(
    log10.(PowerSpectrum.PowerSpectrumNonlinArray),
    BSpline(Cubic(Line(OnGrid()))))
    InterpPmm = scale(InterpPmm, x, y)
    InterpPmm = Interpolations.extrapolate(InterpPmm, Line())
    for iidx in 1:length(ConvolvedDensity.ZBinArray)-1
        for lidx in 1:length(CosmologicalGrid.ℓBinCenters)
            PowerSpectrum.InterpolatedPowerSpectrumBeyondLimber[lidx, :] =
            10 .^(InterpPmm.(log10.(
            CosmologicalGrid.KBeyondLimberArray[lidx, :]),
            CosmologicalGrid.ZArray[1]))
        end
    end
end

function ExtractGrowthFactor!(PowerSpectrum::PowerSpectrum)
    PowerSpectrum.GrowthFactor = zeros(
        length(PowerSpectrum.PowerSpectrumLinArray[1,:]))
    for zidx in 1:length(PowerSpectrum.GrowthFactor)
        PowerSpectrum.GrowthFactor[zidx] =
        sqrt(PowerSpectrum.PowerSpectrumLinArray[1,zidx]./
        PowerSpectrum.PowerSpectrumLinArray[1,1])
    end
end
