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

function CreateDirectories(Cosmologies::Dict, DictCosmo::Dict, path::String)
    mkdir(path)
    mkdir(path*"PowerSpectrum")
    mkdir(path*"Angular")
    mkdir(path*"Derivative/")
    for (key, value) in Cosmologies
        mkdir(path*"Angular/"*key)
        mkdir(path*"PowerSpectrum/"*key)
    end
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
        WriteCosmology(w0waCDMCosmology, PathOutput*key)
    end
end

function InstantiateComputeWeightFunctionOverGrid(
    ConvolvedDensity::AbstractConvolvedDensity,
    w0waCDMCosmology::w0waCDMCosmologyStruct, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    GCWeightFunction::GCWeightFunctionStruct)
    ComputeBiasOverGrid(CosmologicalGrid, GCWeightFunction, ConvolvedDensity)
    ComputeWeightFunctionOverGrid(GCWeightFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
    return GCWeightFunction
end

function InstantiateComputeWeightFunctionOverGrid(
    ConvolvedDensity::AbstractConvolvedDensity,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    WLWeightFunction::WLWeightFunctionStruct)
    ComputeLensingEfficiencyOverGridCustom(WLWeightFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
    ComputeWeightFunctionOverGrid(WLWeightFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
    return WLWeightFunction
end

function EvaluateDerivativeAngularCoefficients(DictCosmo::Dict, Path::String,
    Steps::Array)
    CentralCosmologyCL = ReadAngularCoefficients(
    Path*"/Angular/dvar_central_step_0/cl")
    AngularCoefficientsArray = zeros(size(CentralCosmologyCL.AngularCoefficientsArray, 1),
    size(CentralCosmologyCL.AngularCoefficientsArray, 2),
    size(CentralCosmologyCL.AngularCoefficientsArray, 3),
    2*length(Steps)+1)
    DerivativeArray = similar(CentralCosmologyCL.AngularCoefficientsArray)
    AngularCoefficientsArray[:, :, :, length(Steps) + 1] .=
    CentralCosmologyCL.AngularCoefficientsArray[:,:,:]
    stepvalues = zeros(2*length(Steps)+1)
    AngularDerivatives = zeros(size(AngularCoefficientsArray, 1),
    size(AngularCoefficientsArray, 2), size(AngularCoefficientsArray, 3))
    for (key, value) in DictCosmo
        if value[2] == "present"
            stepvalues[length(Steps) + 1] = value[1]
            for (index, mystep) in enumerate(Steps)
                AngularCoefficientsMinus = ReadAngularCoefficients(
                Path*"/Angular/dvar_"*key*"_step_m_"*string(index)*"/cl")
                AngularCoefficientsPlus = ReadAngularCoefficients(
                Path*"/Angular/dvar_"*key*"_step_p_"*string(index)*"/cl")
                AngularCoefficientsArray[:, :, :, length(Steps) + 1 - index] .=
                AngularCoefficientsMinus.AngularCoefficientsArray
                AngularCoefficientsArray[:, :, :, length(Steps) + 1 + index] .=
                AngularCoefficientsPlus.AngularCoefficientsArray
                stepvalues[length(Steps) + 1 - index] =
                IncrementedValue(value[1], -mystep)
                stepvalues[length(Steps) + 1 + index] =
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

function EvaluateDerivativeAngularCoefficientsNew(DictCosmo::Dict,
    PathInput::String, PathConfig::String, Steps::Array)
    ProbesDict = JSON.parsefile(PathConfig)
    CoefficientsArray = GetProbesArray(ProbesDict)
    for Coefficient in CoefficientsArray
        CentralCosmologyCL = ReadAngularCoefficients(
        PathInput*"/Angular/dvar_central_step_0/cl", Coefficient)
        AngularCoefficientsArray = zeros(size(CentralCosmologyCL.AngularCoefficientsArray, 1),
        size(CentralCosmologyCL.AngularCoefficientsArray, 2),
        size(CentralCosmologyCL.AngularCoefficientsArray, 3),
        2*length(Steps)+1)
        DerivativeArray = similar(CentralCosmologyCL.AngularCoefficientsArray)
        AngularCoefficientsArray[:, :, :, length(Steps) + 1] .=
        CentralCosmologyCL.AngularCoefficientsArray[:,:,:]
        stepvalues = zeros(2*length(Steps)+1)
        AngularDerivatives = zeros(size(AngularCoefficientsArray, 1),
        size(AngularCoefficientsArray, 2), size(AngularCoefficientsArray, 3))
        for (key, value) in DictCosmo
            if value[2] == "present"
                stepvalues[length(Steps) + 1] = value[1]
                for (index, mystep) in enumerate(Steps)
                    AngularCoefficientsMinus = ReadAngularCoefficients(
                    PathInput*"/Angular/dvar_"*key*"_step_m_"*string(index)*"/cl",
                    Coefficient)
                    AngularCoefficientsPlus = ReadAngularCoefficients(
                    PathInput*"/Angular/dvar_"*key*"_step_p_"*string(index)*"/cl",
                    Coefficient)
                    AngularCoefficientsArray[:, :, :, length(Steps) + 1 - index] .=
                    AngularCoefficientsMinus.AngularCoefficientsArray
                    AngularCoefficientsArray[:, :, :, length(Steps) + 1 + index] .=
                    AngularCoefficientsPlus.AngularCoefficientsArray
                    stepvalues[length(Steps) + 1 - index] =
                    IncrementedValue(value[1], -mystep)
                    stepvalues[length(Steps) + 1 + index] =
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
                PathInput*"/Derivative/"*key*"/"*key, Coefficient)
            end
        end
    end
end




function InstantiateWL(DictInput::Dict)
    WLWeightFunction = WLWeightFunctionStruct()
    return WLWeightFunction
end

function InstantiateGC(DictInput::Dict)
    GCWeightFunction = GCWeightFunctionStruct()
    InstantiateBias(DictInput, GCWeightFunction)
    return GCWeightFunction
end

function InstantiateBias(DictInput::Dict,
    GCWeightFunction::GCWeightFunctionStruct)
    if DictInput["PhotometricGalaxy"]["Bias"] == "PiecewiseBias"
        GCWeightFunction.BiasKind = PiecewiseBiasStruct()
    else
        println("Bias must be correctly specified!")
    end
end

function InitializeProbes(DictInput::Dict,
    ConvolvedDensity::AbstractConvolvedDensity,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)
    DictProbes = Dict()
    if DictInput["Lensing"]["present"]
        WLWeightFunction = InstantiateWL(DictInput::Dict)
        WLWeightFunction = InstantiateComputeWeightFunctionOverGrid(ConvolvedDensity,
        w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities,
        WLWeightFunction)
        push!(DictProbes, "Lensing" => WLWeightFunction)
    end
    if DictInput["PhotometricGalaxy"]["present"]
        GCWeightFunction = InstantiateGC(DictInput::Dict)
        GCWeightFunction =
        InstantiateComputeWeightFunctionOverGrid(ConvolvedDensity,
        w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities,
        GCWeightFunction)
        push!(DictProbes, "PhotometricGalaxy" => GCWeightFunction)
    end
    return DictProbes
end

function GetProbesArray(DictInput::Dict)
    ProbesArray = []
    CoefficientsArray = []
    if DictInput["Lensing"]["present"]
        push!(ProbesArray, "Lensing")
    end
    if DictInput["PhotometricGalaxy"]["present"]
        push!(ProbesArray, "PhotometricGalaxy")
    end
    sort!(ProbesArray)
    for key_A in ProbesArray
        for key_B in ProbesArray
            if key_B*"_"*key_A in CoefficientsArray
            else
                push!(CoefficientsArray, key_A*"_"*key_B)
            end
        end
    end
    return CoefficientsArray
end

function InitializeComputeAngularCoefficients(ProbesDict::Dict,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    CosmologicalGrid::CosmologicalGrid, PowerSpectrum::PowerSpectrum,
    PathOutput::String, key::String)
    ProbesArray = []
    CoefficientsArray = []
    for (key, value) in ProbesDict
        push!(ProbesArray, key)
    end
    sort!(ProbesArray)
    for key_A in ProbesArray
        for key_B in ProbesArray
            if key_B*"_"*key_A in CoefficientsArray
            else
                push!(CoefficientsArray, key_A*"_"*key_B)
                AngularCoefficients = AngularCoefficientsStruct(
                AngularCoefficientsArray
                = zeros(length(CosmologicalGrid.MultipolesArray),
                length(ProbesDict[key_A].WeightFunctionArray[:, 1]),
                length(ProbesDict[key_B].WeightFunctionArray[:, 1])))
                ComputeAngularCoefficients(AngularCoefficients,
                ProbesDict[key_A], ProbesDict[key_B], BackgroundQuantities,
                w0waCDMCosmology, CosmologicalGrid, PowerSpectrum,
                CosmoCentral.CustomTrapz())
                WriteAngularCoefficients(key_A*"_"*key_B,
                AngularCoefficients, PathOutput*key*"/cl")
                WriteCosmology(w0waCDMCosmology, PathOutput*key)
            end
        end
    end
end

function InitializeComputeAngularCoefficients(ProbesDict::Dict,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    CosmologicalGrid::CosmologicalGrid, PowerSpectrum::PowerSpectrum,
    PathOutput::String, CosmoDict::Dict, key::String)
    ProbesArray = []
    CoefficientsArray = []
    for (key, value) in ProbesDict
        push!(ProbesArray, key)
    end
    sort!(ProbesArray)
    for key_A in ProbesArray
        for key_B in ProbesArray
            if key_B*"_"*key_A in CoefficientsArray
            else
                push!(CoefficientsArray, key_A*"_"*key_B)
                AngularCoefficients = AngularCoefficientsStruct(
                AngularCoefficientsArray
                = zeros(length(CosmologicalGrid.MultipolesArray),
                length(ProbesDict[key_A].WeightFunctionArray[:, 1]),
                length(ProbesDict[key_B].WeightFunctionArray[:, 1])))
                ComputeAngularCoefficients(AngularCoefficients,
                ProbesDict[key_A], ProbesDict[key_B], BackgroundQuantities,
                w0waCDMCosmology, CosmologicalGrid, PowerSpectrum,
                CosmoCentral.CustomTrapz())
                WriteAngularCoefficients(key_A*"_"*key_B,
                AngularCoefficients, PathOutput*key*"/cl")
                WriteParameters(CosmoDict, PathOutput*key)
            end
        end
    end
end

function EvaluateAngularCoefficients(Cosmologies::Dict, PathInput::String,
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
