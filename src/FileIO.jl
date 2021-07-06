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

function ReadCosmologyForecast(CosmoDict::Dict)
    Cosmology = w0waCDMCosmology(
    w0 = CosmoDict["w0"][1],
    wa = CosmoDict["wa"][1],
    Mν = CosmoDict["Mν"][1],
    H0 = CosmoDict["H0"][1],
    ΩM = CosmoDict["ΩM"][1],
    ΩB = CosmoDict["ΩB"][1],
    ΩDE = CosmoDict["ΩDE"][1],
    Ωk = CosmoDict["Ωk"][1],
    Ωr = CosmoDict["Ωr"][1],
    ns = CosmoDict["ns"][1],
    σ8 = CosmoDict["σ8"][1])
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

function ReadIntrinsicAlignment(CosmoDict::Dict)
    intrinsicalignment = ExtendedNLIA()
    intrinsicalignment.𝓐IA = CosmoDict["𝓐IA"]
    intrinsicalignment.βIA = CosmoDict["βIA"]
    intrinsicalignment.𝓒IA = CosmoDict["𝓒IA"]
    intrinsicalignment.ηIA = CosmoDict["ηIA"]
    return intrinsicalignment
end


function ReadIntrinsicAlignmentForecast(IADict::Dict, IAModel::String)
    if IAModel == "ExtendedNLIA"
        intrinsicalignment = ExtendedNLIA()
        intrinsicalignment.𝓐IA = IADict["𝓐IA"][1]
        intrinsicalignment.βIA = IADict["βIA"][1]
        intrinsicalignment.𝓒IA = IADict["𝓒IA"][1]
        intrinsicalignment.ηIA = IADict["ηIA"][1]
    elseif IAModel == "None"
        intrinsicalignment = AbsentIA()
    else
        ErrorException("Intrinsc Alignment model not defined correctly.")
    end
    return intrinsicalignment
end

function ReadBias(CosmoDict::Dict)
    bias = EuclidBias()
    bias.A = CosmoDict["A"]
    bias.B = CosmoDict["B"]
    bias.C = CosmoDict["C"]
    bias.D = CosmoDict["D"]
    return bias
end

function ReadBiasForecast(BiasDict::Dict, BiasModel::String)
    if BiasModel == "EuclidBias"
        bias = EuclidBias()
        bias.A = BiasDict["A"][1]
        bias.B = BiasDict["B"][1]
        bias.C = BiasDict["C"][1]
        bias.D = BiasDict["D"][1]
    elseif BiasModel == "PiecewiseBias"
        bias = PiecewiseBias()
        bias.BiasMultiplier = BiasDict["BiasMultiplier"]
    else
        ErrorException("Bias model not defined correctly.")
    end
    return bias
end

function ReadCℓ(Filename::String, Probes::String)
    Filename *= ".h5"
    file = HDF5.h5open(Filename, "r")
    c_lij = HDF5.read(file["cls"][Probes]["c_lij"])
    return Cℓ(CℓArray = c_lij)
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
