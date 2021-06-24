function WriteCℓ!(Probes::String,
    Cℓ::AbstractCℓ, Filename::String)
    h5write(Filename*".h5", "cls/"*Probes*"/c_lij",
    Cℓ.CℓArray)
end

function WriteCosmology!(Cosmology::w0waCDMCosmology, Filename::String)
    CosmoDict = Dict{String,Float64}()
    CosmoDict["w0"]  = Cosmology.w0
    CosmoDict["wa"]  = Cosmology.wa
    CosmoDict["Mν"]  = Cosmology.Mν
    CosmoDict["H0"]  = Cosmology.H0
    CosmoDict["ΩM"]  = Cosmology.ΩM
    CosmoDict["ΩB"]  = Cosmology.ΩB
    CosmoDict["ΩDE"] = Cosmology.ΩDE
    CosmoDict["Ωk"]  = Cosmology.Ωk
    CosmoDict["Ωr"]  = Cosmology.Ωr
    CosmoDict["ns"]  = Cosmology.ns
    CosmoDict["σ8"]  = Cosmology.σ8
    JSON3.write(Filename*"/cosmology.json", CosmoDict)
end

function WriteParameters!(CosmoDict::Dict, Filename::String)
    JSON3.write(Filename*"/parameters.json", CosmoDict)
end

function ReadCosmology(CosmoDict::Dict)
    Cosmology = w0waCDMCosmology(
    w0 = CosmoDict["w0"],
    wa = CosmoDict["wa"],
    Mν = CosmoDict["Mν"],
    H0 = CosmoDict["H0"],
    ΩM = CosmoDict["ΩM"],
    ΩB = CosmoDict["ΩB"],
    ΩDE = CosmoDict["ΩDE"],
    Ωk = CosmoDict["Ωk"],
    Ωr = CosmoDict["Ωr"],
    ns = CosmoDict["ns"],
    σ8 = CosmoDict["σ8"])
    return Cosmology
end

function ReadCosmology(CosmoDict::JSON3.Object)
    Cosmology = w0waCDMCosmology(
    w0 = CosmoDict["w0"],
    wa = CosmoDict["wa"],
    Mν = CosmoDict["Mν"],
    H0 = CosmoDict["H0"],
    ΩM = CosmoDict["ΩM"],
    ΩB = CosmoDict["ΩB"],
    ΩDE = CosmoDict["ΩDE"],
    Ωk = CosmoDict["Ωk"],
    Ωr = CosmoDict["Ωr"],
    ns = CosmoDict["ns"],
    σ8 = CosmoDict["σ8"])
    return Cosmology
end

function ReadCℓ(Filename::String)
    Filename *= ".h5"
    file = HDF5.h5open(Filename, "r")
    c_lij =
    HDF5.read(file["cls"]["PhotometricGalaxy_PhotometricGalaxy"]["c_lij"])
    AngularCoefficients = Cℓ(CℓArray =
    c_lij)
    return AngularCoefficients
end

function ReadCℓ(Filename::String, Probes::String)
    Filename *= ".h5"
    file = HDF5.h5open(Filename, "r")
    c_lij =
    HDF5.read(file["cls"][Probes]["c_lij"])
    AngularCoefficients = AngularCoefficients(CℓArray =
    c_lij)
    return AngularCoefficients
end

function Write∂Cℓ!(DerivativeArray::AbstractArray{Float64, 3},
    Filename::String)
    h5write(Filename*".h5", "dcls/PhotometricGalaxy_PhotometricGalaxy/dc_lij",
    DerivativeArray)
end

function Write∂Cℓ!(DerivativeArray::AbstractArray{Float64, 3},
    Filename::String, Probes::String)
    h5write(Filename*".h5", "dcls/"*Probes*"/dc_lij",
    DerivativeArray)
end


"""
    WritePowerSpectrumBackground(PowerSpectrum::PowerSpectrum,
    BackgroundQuantities::BackgroundQuantities,
    CosmologicalGrid::CosmologicalGrid, Filename::String)

This function writes the Power Spectrum, the Background quantities and the
Cosmological Grid in a HDF5 file.
"""
function WritePowerSpectrumBackground(PowerSpectrum::AbstractPowerSpectrum,
    BackgroundQuantities::AbstractBackgroundQuantities,
    CosmologicalGrid::AbstractCosmologicalGrid, Filename::String)
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

"""
    ReadPowerSpectrumBackground(Filename::String,
    MultipolesArray::Vector{Float64})

This function reads the Power Spectrum, the Background quantities and the
Cosmological Grid from a HDF5 file.
"""
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
    BackgroundQuantitiesRead = BackgroundQuantities(HZArray = hubble_array,
    rZArray = comoving_distance_array)
    CosmologicalGridRead = CosmologicalGrid(ZArray = z_grid, KArray = k_grid,
    MultipolesArray = MultipolesArray)
    PowerSpectrumRead = PowerSpectrum(PowerSpectrumLinArray = lin_p_mm_k_z,
    PowerSpectrumNonlinArray = nonlin_p_mm_k_z,
    InterpolatedPowerSpectrum = zeros(length(CosmologicalGridRead.MultipolesArray),
    length(CosmologicalGridRead.ZArray)))
    return PowerSpectrumRead, BackgroundQuantitiesRead, CosmologicalGridRead
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
    BackgroundQuantities = BackgroundQuantities(HZArray =
    dimensionless_hubble_array*67,
    rZArray = dimensionless_comoving_distance_array*c_0/67)
    CosmologicalGrid = CosmologicalGrid(ZArray = z_grid, KArray = k_grid,
    MultipolesArray = MultipolesArray)
    PowerSpectrum = PowerSpectrum(PowerSpectrumLinArray = lin_p_mm_k_z,
    PowerSpectrumNonlinArray = nonlin_p_mm_k_z,
    InterpolatedPowerSpectrum = zeros(length(CosmologicalGrid.MultipolesArray),
    length(CosmologicalGrid.ZArray)))
    return PowerSpectrum, BackgroundQuantities, CosmologicalGrid
end
