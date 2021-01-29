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
                myvalue = IncrementedValue(value[1], -mystep)
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

function CreateDirectoriesDerivatives(DictCosmo::Dict, path::String)
    mkdir(path*"Derivative/")
    for (key, value) in DictCosmo
        if value[2] == "present"
            mkdir(path*"Derivative/"*key)
        end
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
        PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
        ReadPowerSpectrumBackground(PathInput*key*"/p_mm",
        CosmologicalGrid.MultipolesArray)
        GCWeightFunction = GCWeightFunctionStruct()
        WLWeightFunction = WLWeightFunctionStruct()
        GCWeightFunction =
        InstantiateComputeWeightFunctionOverGrid(ConvolvedDensity,
        w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities,
        GCWeightFunction)
        ComputeLimberArray(CosmologicalGrid, BackgroundQuantities)
        WLWeightFunction =
        InstantiateComputeWeightFunctionOverGrid(ConvolvedDensity,
        w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities,
        WLWeightFunction)
        InterpolateAndEvaluatePowerSpectrum(CosmologicalGrid,
        BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
        GCGCAngularCoefficients = GCGCAngularCoefficientsStruct(
        AngularCoefficientsArray
        = zeros(length(CosmologicalGrid.MultipolesArray),
        length(GCWeightFunction.WeightFunctionArray[:, 1]),
        length(GCWeightFunction.WeightFunctionArray[:, 1])))
        ComputeAngularCoefficients(GCGCAngularCoefficients, GCWeightFunction,
        GCWeightFunction, BackgroundQuantities, w0waCDMCosmology,
        CosmologicalGrid, PowerSpectrum,
        CosmoCentral.CustomTrapz())
        WriteAngularCoefficients(GCGCAngularCoefficients, CosmologicalGrid,
        GCWeightFunction, ConvolvedDensity, PathOutput*key*"/cl")
        WLWLAngularCoefficients = WLWLAngularCoefficientsStruct(AngularCoefficientsArray
        = zeros(length(CosmologicalGrid.MultipolesArray),
        length(WLWeightFunction.WeightFunctionArray[:, 1]),
        length(WLWeightFunction.WeightFunctionArray[:, 1])))
        ComputeAngularCoefficients(WLWLAngularCoefficients, WLWeightFunction,
        WLWeightFunction, BackgroundQuantities, w0waCDMCosmology,
        CosmologicalGrid, PowerSpectrum,
        CosmoCentral.CustomTrapz())
        WriteAngularCoefficients(WLWLAngularCoefficients, CosmologicalGrid,
        WLWeightFunction, ConvolvedDensity, PathOutput*key*"/cl")
        GCWLAngularCoefficients = GCWLAngularCoefficientsStruct(AngularCoefficientsArray
        = zeros(length(CosmologicalGrid.MultipolesArray),
        length(WLWeightFunction.WeightFunctionArray[:, 1]),
        length(WLWeightFunction.WeightFunctionArray[:, 1])))
        ComputeAngularCoefficients(GCWLAngularCoefficients, GCWeightFunction,
        WLWeightFunction, BackgroundQuantities, w0waCDMCosmology,
        CosmologicalGrid, PowerSpectrum,
        CosmoCentral.CustomTrapz())
        WriteAngularCoefficients(GCWLAngularCoefficients, CosmologicalGrid,
        ConvolvedDensity, PathOutput*key*"/cl")
    end
end

function InstantiateComputeWeightFunctionOverGrid(
    ConvolvedDensity::AsbtractConvolvedDensity,
    w0waCDMCosmology::w0waCDMCosmologyStruct, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    GCWeightFunction::GCWeightFunctionStruct)
    #PiecewiseBias = PiecewiseBiasStruct(BiasArray =
    #zeros(length(ConvolvedDensity.DensityNormalizationArray),
    #length(CosmologicalGrid.ZArray)))
    GCWeightFunction = GCWeightFunctionStruct(WeightFunctionArray =
    zeros(length(ConvolvedDensity.DensityNormalizationArray),
    length(CosmologicalGrid.ZArray)))
    ComputeBiasOverGrid(CosmologicalGrid, GCWeightFunction,
    GCWeightFunction.BiasKind,
    ConvolvedDensity)
    ComputeWeightFunctionOverGrid(GCWeightFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
    return GCWeightFunction
end

function InstantiateComputeWeightFunctionOverGrid(ConvolvedDensity::AsbtractConvolvedDensity,
    w0waCDMCosmology::w0waCDMCosmologyStruct, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    WLWeightFunction::WLWeightFunctionStruct)
    WLWeightFunction = WLWeightFunctionStruct(WeightFunctionArray =
    zeros(length(ConvolvedDensity.DensityNormalizationArray),
    length(CosmologicalGrid.ZArray)), LensingEfficiencyArray =
    zeros(length(ConvolvedDensity.DensityNormalizationArray),
    length(CosmologicalGrid.ZArray)))
    ComputeLensingEfficiencyOverGridCustom(WLWeightFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
    ComputeWeightFunctionOverGrid(WLWeightFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
    return WLWeightFunction
end

function EvaluateDerivativeAngularCoefficients(DictCosmo::Dict, Path::String,
    steps::Array)
    CentralCosmologyCL = ReadAngularCoefficients(
    Path*"/Angular/dvar_central_step_0/cl")
    AngularCoefficientsArray = zeros(size(CentralCosmologyCL.AngularCoefficientsArray, 1),
    size(CentralCosmologyCL.AngularCoefficientsArray, 2),
    size(CentralCosmologyCL.AngularCoefficientsArray, 3),
    2*length(steps)+1)
    DerivativeArray = similar(CentralCosmologyCL.AngularCoefficientsArray)
    AngularCoefficientsArray[:, :, :, length(steps) + 1] .=
    CentralCosmologyCL.AngularCoefficientsArray[:,:,:]
    stepvalues = zeros(2*length(steps)+1)
    AngularDerivatives = zeros(size(AngularCoefficientsArray, 1),
    size(AngularCoefficientsArray, 2), size(AngularCoefficientsArray, 3))
    for (key, value) in DictCosmo
        if value[2] == "present"
            stepvalues[length(steps) + 1] = value[1]
            for (index, mystep) in enumerate(steps)
                AngularCoefficientsMinus = ReadAngularCoefficients(
                Path*"/Angular/dvar_"*key*"_step_m_"*string(index)*"/cl")
                AngularCoefficientsPlus = ReadAngularCoefficients(
                Path*"/Angular/dvar_"*key*"_step_p_"*string(index)*"/cl")
                AngularCoefficientsArray[:, :, :, length(steps) + 1 - index] .=
                AngularCoefficientsMinus.AngularCoefficientsArray
                AngularCoefficientsArray[:, :, :, length(steps) + 1 + index] .=
                AngularCoefficientsPlus.AngularCoefficientsArray
                stepvalues[length(steps) + 1 - index] =
                IncrementedValue(value[1], -mystep)
                stepvalues[length(steps) + 1 + index] =
                IncrementedValue(value[1],  mystep)
            end
            for idx_a in 1:size(CentralCosmologyCL.AngularCoefficientsArray, 2)
                for idx_b in 1:size(CentralCosmologyCL.AngularCoefficientsArray, 3)
                    for idx_l in 1:size(CentralCosmologyCL.AngularCoefficientsArray, 1)
                        y = AngularCoefficientsArray[idx_l, idx_a, idx_b, :]
                        der = SteMDerivative(stepvalues, y)
                        AngularDerivatives[idx_l, idx_a, idx_b] = der
                    end
                end
            end
            WriteDerivativeCoefficients(AngularDerivatives,
            Path*"/Derivative/"*key*"/"*key)
        end
    end
end
