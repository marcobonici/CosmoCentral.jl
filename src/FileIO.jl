function WriteAngularCoefficients(AngularCoefficients::AngularCoefficients,
    CosmologicalGrid::CosmologicalGrid, WeightFunction::WeightFunction,
    Bias::Bias, ConvolvedDensity::ConvolvedDensity, Filename::String)
    h5write(Filename*".h5", "cls/PhotometricGalaxy_PhotometricGalaxy/c_lij",
    AngularCoefficients.AngularCoefficientsArray)
    h5write(Filename*".h5",
    "weight_functions/PhotometricGalaxy/bias/b_i_z",
    Bias.BiasArray)
    h5write(Filename*".h5",
    "weight_functions/PhotometricGalaxy/density/norm_density_iz",
    ConvolvedDensity.DensityNormalizationArray)
end

function WritePowerSpectrumBackground(PowerSpectrum::PowerSpectrum,
    BackgroundQuantities::BackgroundQuantities,
    CosmologicalGrid::CosmologicalGrid, Filename::String)
    h5write(Filename*".h5",
    "cosmology/comoving_distance_array",
    BackgroundQuantities.rZArray)
    h5write(Filename*".h5",
    "cosmology/hubble_array",
    BackgroundQuantities.HZArray)
    h5write(Filename*".h5",
    "cosmology/z_grid",
    CosmologicalGrid.ZArray)
    h5write(Filename*".h5",
    "power_spectrum/z_grid",
    CosmologicalGrid.ZArray)
    h5write(Filename*".h5",
    "power_spectrum/k_grid",
    CosmologicalGrid.KArray)
    h5write(Filename*".h5",
    "power_spectrum/lin_p_mm_z_k",
    PowerSpectrum.PowerSpectrumLinArray)
    h5write(Filename*".h5",
    "power_spectrum/nonlin_p_mm_z_k",
    PowerSpectrum.PowerSpectrumNonlinArray)
end
