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
