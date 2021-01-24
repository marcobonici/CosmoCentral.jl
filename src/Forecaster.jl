function IncrementedValue(value, step)
    increment = value
    if increment == 0
        increment = 1.
    end
    myvalue = value + abs(increment) * (step)
    return myvalue
end

function CreateCosmologies(DictCosmo::Dict, steps::Array)
    MyDict = Dict{String,Array{Any,1}}()
    for (key, value) in DictCosmo
        if value[2] == "present"
            for (index, mystep) in enumerate(steps)
                CopyDictCosmo = deepcopy(DictCosmo)
                myvalue = IncrementedValue(value[1], mystep)
                CopyDictCosmo[key] = [myvalue]
                w0waCDMCosmology = CosmoCentral.w0waCDMCosmologyStruct(
                w0 = CopyDictCosmo["w0"][1], wa = CopyDictCosmo["wa"][1],
                Mν = CopyDictCosmo["Mν"][1], H0 = CopyDictCosmo["H0"][1],
                ΩM = CopyDictCosmo["ΩM"][1], ΩB = CopyDictCosmo["ΩB"][1],
                ΩDE = CopyDictCosmo["ΩDE"][1], Ωk = CopyDictCosmo["Ωk"][1],
                Ωr = CopyDictCosmo["Ωr"][1], ns = CopyDictCosmo["ns"][1],
                σ8 = CopyDictCosmo["σ8"][1])
                MyDict["dvar_"*key*"_step_p_"*string(index)] = [w0waCDMCosmology]
                myvalue = IncrementedValue(value[1], mystep)
                CopyDictCosmo[key] = [myvalue]
                w0waCDMCosmology = CosmoCentral.w0waCDMCosmologyStruct(w0 =
                CopyDictCosmo["w0"][1], wa = CopyDictCosmo["wa"][1],
                Mν = CopyDictCosmo["Mν"][1], H0 = CopyDictCosmo["H0"][1],
                ΩM = CopyDictCosmo["ΩM"][1], ΩB = CopyDictCosmo["ΩB"][1],
                ΩDE = CopyDictCosmo["ΩDE"][1], Ωk = CopyDictCosmo["Ωk"][1],
                Ωr = CopyDictCosmo["Ωr"][1], ns = CopyDictCosmo["ns"][1],
                σ8 = CopyDictCosmo["σ8"][1])
                MyDict["dvar_"*key*"_step_m_"*string(index)] = [w0waCDMCosmology]
            end
        end
    end
    w0waCDMCosmology = CosmoCentral.w0waCDMCosmologyStruct(
    w0 = DictCosmo["w0"][1], wa = DictCosmo["wa"][1], Mν = DictCosmo["Mν"][1],
    H0 = DictCosmo["H0"][1], ΩM = DictCosmo["ΩM"][1], ΩB = DictCosmo["ΩB"][1],
    ΩDE = DictCosmo["ΩDE"][1], Ωk = DictCosmo["Ωk"][1], Ωr = DictCosmo["Ωr"][1],
    ns = DictCosmo["ns"][1], σ8 = DictCosmo["σ8"][1])
    MyDict["dvar_central_step_0"] = [w0waCDMCosmology]
    return MyDict
end

function CreateDirectories(Cosmologies::Dict, path::String)
    mkdir(path)
    mkdir(path*"PowerSpectrum")
    mkdir(path*"Angular")
    for (key, value) in Cosmologies
    mkdir(path*"Angular/"*key)
    mkdir(path*"PowerSpectrum/"*key)
    end
end

function EvaluatePowerSpectra(Cosmologies::Dict, Path::String,
    CosmologicalGrid::CosmologicalGrid)
    for (key, value) in Cosmologies
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
    end
end

function EvaluateAngularCoefficients(Cosmologies::Dict, PathInput::String,
    PathOutput::String, CosmologicalGrid::CosmologicalGrid)
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
        MultipolesArray = Array(LogSpaced(10., 3000., 100))
        PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
        ReadPowerSpectrumBackground(PathInput*key*"/p_mm", MultipolesArray)
        GCWeightFunction, PiecewiseBias =
        InstantiateComputeWeightFunctionOverGrid(AnalitycalDensity,
        InstrumentResponse, ConvolvedDensity, w0waCDMCosmology,
        CosmologicalGrid, BackgroundQuantities)
        ComputeLimberArray(CosmologicalGrid, BackgroundQuantities)
        InterpolateAndEvaluatePowerSpectrum(CosmologicalGrid,
        BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
        AngularCoefficients = AngularCoefficientsStruct(AngularCoefficientsArray
        = zeros(length(CosmologicalGrid.MultipolesArray),
        length(GCWeightFunction.WeightFunctionArray[:, 1]),
        length(GCWeightFunction.WeightFunctionArray[:, 1])))
        ComputeAngularCoefficients(AngularCoefficients, GCWeightFunction,
        GCWeightFunction, BackgroundQuantities, w0waCDMCosmology,
        CosmologicalGrid, PowerSpectrum,
        CosmoCentral.NumericalIntegrationSimpson())
        WriteAngularCoefficients(AngularCoefficients, CosmologicalGrid,
        GCWeightFunction, PiecewiseBias, ConvolvedDensity, PathOutput*key*"/cl")
    end
end

function InstantiateComputeWeightFunctionOverGrid(
    AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, ConvolvedDensity::ConvolvedDensity,
    w0waCDMCosmology::w0waCDMCosmology, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)
    PiecewiseBias = PiecewiseBiasStruct(BiasArray =
    zeros(length(ConvolvedDensity.DensityNormalizationArray),
    length(CosmologicalGrid.ZArray)))
    ComputeBiasOverGrid(CosmologicalGrid, PiecewiseBias, ConvolvedDensity)
    GCWeightFunction = GCWeightFunctionStruct(WeightFunctionArray =
    zeros(length(ConvolvedDensity.DensityNormalizationArray),
    length(CosmologicalGrid.ZArray)))
    ComputeWeightFunctionOverGrid(GCWeightFunction, AnalitycalDensity,
    InstrumentResponse, ConvolvedDensity, PiecewiseBias, CosmologicalGrid,
    BackgroundQuantities, w0waCDMCosmology)
    return GCWeightFunction, PiecewiseBias
end
