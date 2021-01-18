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

function ReadAngularCoefficients(Filename::String)
    Filename *= ".h5"
    file = HDF5.h5open(Filename, "r")
    c_lij =
    HDF5.read(file["cls"]["PhotometricGalaxy_PhotometricGalaxy"]["c_lij"])
    AngularCoefficients = AngularCoefficientsStruct(AngularCoefficientsArray =
    c_lij)
    return AngularCoefficients
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
    "power_spectrum/lin_p_mm_k_z",
    PowerSpectrum.PowerSpectrumLinArray)
    h5write(Filename*".h5",
    "power_spectrum/nonlin_p_mm_k_z",
    PowerSpectrum.PowerSpectrumNonlinArray)
end

function ReadPowerSpectrumBackground(Filename::String,
    MultipolesArray::Vector{Float64})
    Filename *= ".h5"
    file = HDF5.h5open(Filename, "r")
    nonlin_p_mm_k_z = HDF5.read(file["power_spectrum"]["nonlin_p_mm_k_z"])
    lin_p_mm_k_z    = HDF5.read(file["power_spectrum"]["lin_p_mm_k_z"])
    z_grid = HDF5.read(file["power_spectrum"]["z_grid"])
    k_grid = HDF5.read(file["power_spectrum"]["k_grid"])
    comoving_distance_array = HDF5.read(file["cosmology"]["comoving_distance_array"])
    hubble_array = HDF5.read(file["cosmology"]["hubble_array"])
    BackgroundQuantities = BackgroundQuantitiesStruct(HZArray = hubble_array,
    rZArray = comoving_distance_array)
    CosmologicalGrid = CosmologicalGridStruct(ZArray = z_grid, KArray = k_grid,
    MultipolesArray = MultipolesArray)
    PowerSpectrum = PowerSpectrumStruct(PowerSpectrumLinArray = lin_p_mm_k_z,
    PowerSpectrumNonlinArray = nonlin_p_mm_k_z,
    InterpolatedPowerSpectrum = zeros(length(CosmologicalGrid.MultipolesArray), length(CosmologicalGrid.ZArray)))
    return PowerSpectrum, BackgroundQuantities, CosmologicalGrid
end

function ReadPowerSpectrumBackgroundSeyfert(Filename::String,
    MultipolesArray::Vector{Float64})
    c_0 = 2.99792458e5
    Filename *= ".h5"
    file = HDF5.h5open(Filename, "r")
    nonlin_p_mm_k_z = HDF5.read(file["power_spectrum"]["nonlin_p_mm_z_k"])
    lin_p_mm_k_z    = HDF5.read(file["power_spectrum"]["lin_p_mm_z_k"])
    z_grid = HDF5.read(file["power_spectrum"]["z_grid"])
    k_grid = HDF5.read(file["power_spectrum"]["k_grid"])
    dimensionless_comoving_distance_array = HDF5.read(file["cosmology"]["dimensionless_comoving_distance_array"])
    dimensionless_hubble_array = HDF5.read(file["cosmology"]["dimensionless_hubble_array"])
    BackgroundQuantities = BackgroundQuantitiesStruct(HZArray =
    dimensionless_hubble_array*67,
    rZArray = dimensionless_comoving_distance_array*c_0/67)
    CosmologicalGrid = CosmologicalGridStruct(ZArray = z_grid, KArray = k_grid,
    MultipolesArray = MultipolesArray)
    PowerSpectrum = PowerSpectrumStruct(PowerSpectrumLinArray = lin_p_mm_k_z,
    PowerSpectrumNonlinArray = nonlin_p_mm_k_z,
    InterpolatedPowerSpectrum = zeros(length(CosmologicalGrid.MultipolesArray), length(CosmologicalGrid.ZArray)))
    return PowerSpectrum, BackgroundQuantities, CosmologicalGrid
end

function ReadAngularCoefficients(Filename::String,
    MultipolesArray::Vector{Float64})
    Filename *= ".h5"
    file = HDF5.h5open(Filename, "r")
    nonlin_p_mm_k_z = HDF5.read(file["power_spectrum"]["nonlin_p_mm_k_z"])
    lin_p_mm_k_z    = HDF5.read(file["power_spectrum"]["lin_p_mm_k_z"])
    z_grid = HDF5.read(file["power_spectrum"]["z_grid"])
    k_grid = HDF5.read(file["power_spectrum"]["k_grid"])
    comoving_distance_array = HDF5.read(file["cosmology"]["comoving_distance_array"])
    hubble_array = HDF5.read(file["cosmology"]["hubble_array"])
    BackgroundQuantities = BackgroundQuantitiesStruct(HZArray = hubble_array,
    rZArray = comoving_distance_array)
    CosmologicalGrid = CosmologicalGridStruct(ZArray = z_grid, KArray = k_grid,
    MultipolesArray = MultipolesArray)
    PowerSpectrum = PowerSpectrumStruct(PowerSpectrumLinArray = lin_p_mm_k_z,
    PowerSpectrumNonlinArray = nonlin_p_mm_k_z,
    InterpolatedPowerSpectrum = zeros(length(CosmologicalGrid.MultipolesArray), length(CosmologicalGrid.ZArray)))
    return PowerSpectrum, BackgroundQuantities, CosmologicalGrid
end
