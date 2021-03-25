function CosmologicalParameterSampler(DictCosmo::Dict, n_points::Int)
    parameter_list = ["w0", "wa", "Mν", "H0","ΩM","ΩB","ΩDE","Ωk","Ωr","ns",
    "σ8"]
    dim_space = 0
    priors = [(-1000.0,1.0)]
    deleteat!(priors, 1)
    for key in parameter_list
        if DictCosmo[key][3] == "present"
            dim_space += 1
            append!(priors, Tuple[(DictCosmo[key][1],DictCosmo[key][2])])
        end
    end
    plan, _ = LatinHypercubeSampling.LHCoptim(n_points, dim_space, 2000)
    scaled_plan = LatinHypercubeSampling.scaleLHC(plan, priors)
    CosmoStructDict = Dict{String,Array{Any,1}}()
    for i in 1:n_points
        CopyDictCosmo = deepcopy(DictCosmo)
        counter = 0
        for key in parameter_list
            if DictCosmo[key][3] == "present"
                counter += 1
                CopyDictCosmo[key][4] = scaled_plan[i, counter]
            end
        end
        w0waCDMCosmology = CosmoCentral.w0waCDMCosmologyStruct(
        w0 = CopyDictCosmo["w0"][4], wa = CopyDictCosmo["wa"][4],
        Mν = CopyDictCosmo["Mν"][4], H0 = CopyDictCosmo["H0"][4],
        ΩM = CopyDictCosmo["ΩM"][4], ΩB = CopyDictCosmo["ΩB"][4],
        ΩDE = 1 - CopyDictCosmo["ΩM"][4], Ωk = CopyDictCosmo["Ωk"][4],
        Ωr = CopyDictCosmo["Ωr"][4], ns = CopyDictCosmo["ns"][4],
        σ8 = CopyDictCosmo["σ8"][4])
        CosmoStructDict[string(i)] = [w0waCDMCosmology]
    end
    return scaled_plan, CosmoStructDict
end

function EvaluatePowerSpectraLHS(Cosmologies::Dict, Path::String,
    CosmologicalGrid::CosmologicalGrid)
    for (key, value) in Cosmologies
        println(key, value[1])
        BackgroundQuantities = BackgroundQuantitiesStruct(
        HZArray = zeros(length(CosmologicalGrid.ZArray)),
        rZArray=zeros(length(CosmologicalGrid.ZArray)))
        ComputeBackgroundQuantitiesOverGrid(CosmologicalGrid,
        BackgroundQuantities, value[1])
        ClassyParams = Initializeclassy(value[1])
        PowerSpectrum = PowerSpectrumStruct(PowerSpectrumLinArray =
        zeros(length(CosmologicalGrid.KArray), length(CosmologicalGrid.ZArray)),
        PowerSpectrumNonlinArray = zeros(length(CosmologicalGrid.KArray),
        length(CosmologicalGrid.ZArray)),
        InterpolatedPowerSpectrum = zeros(length(
        CosmologicalGrid.MultipolesArray), length(CosmologicalGrid.ZArray)))
        EvaluatePowerSpectrum(ClassyParams, CosmologicalGrid, PowerSpectrum)
        WritePowerSpectrumBackground(PowerSpectrum, BackgroundQuantities,
        CosmologicalGrid, Path*key*"/p_mm")
        WriteCosmology(value[1], Path*key)
    end
end

function CreateDirectoriesLHS(Cosmologies::Dict, path::String)
    mkdir(path)
    mkdir(path*"PowerSpectrum")
    mkdir(path*"Angular")
    for (key, value) in Cosmologies
        mkdir(path*"Angular/"*key)
        mkdir(path*"PowerSpectrum/"*key)
    end
end

function EvaluateAngularCoefficientsLHS(Cosmologies::Dict, PathInput::String,
    PathOutput::String, CosmologicalGrid::CosmologicalGrid, PathConfig::String)
    ProbesDict = JSON.parsefile(PathConfig)
    AnalitycalDensity = AnalitycalDensityStruct()
    NormalizeAnalitycalDensityStruct(AnalitycalDensity)
    InstrumentResponse = InstrumentResponseStruct()
    ConvolvedDensity = ConvolvedDensityStruct(DensityGridArray =
    ones(10, length(CosmologicalGrid.ZArray)))
    NormalizeConvolvedDensityStruct(ConvolvedDensity, AnalitycalDensity,
    InstrumentResponse, CosmologicalGrid)
    ComputeConvolvedDensityFunctionGrid(CosmologicalGrid, ConvolvedDensity,
    AnalitycalDensity, InstrumentResponse)
    for (key, value) in Cosmologies
        w0waCDMCosmology = value[1]
        PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
        ReadPowerSpectrumBackground(PathInput*key*"/p_mm",
        CosmologicalGrid.MultipolesArray)
        CosmologicalGrid.MultipolesArray[1,1]
        DictProbes = InitializeProbes(ProbesDict, ConvolvedDensity,
        w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities)
        ComputeLimberArray(CosmologicalGrid, BackgroundQuantities)
        InterpolateAndEvaluatePowerSpectrum(CosmologicalGrid,
        BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
        InitializeComputeAngularCoefficients(DictProbes, BackgroundQuantities,
        w0waCDMCosmology, CosmologicalGrid, PowerSpectrum, PathOutput, key)
    end
end