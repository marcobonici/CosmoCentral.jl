function IncrementedValue(value, step)
    increment = value
    if increment == 0
        increment = 1.
    end
    myvalue = value + abs(increment) * (step)
    return myvalue
end


macro IterateParameters(ReferenceDict, VariedDict, ReadMethod, Steps, Model)
    quote
        for (key, value) in $(esc(ReferenceDict))
            if value[2] == "present"
                for (index, mystep) in enumerate($(esc(Steps)))
                    CopyDictCosmo = deepcopy($(esc(ReferenceDict)))
                    myvalue = IncrementedValue(value[1], mystep)
                    CopyDictCosmo[key] = [myvalue]
                    iterated = $(ReadMethod)(CopyDictCosmo, $(esc(Model)))
                    $(esc(VariedDict))["dvar_"*key*"_step_p_"*string(index)] = [iterated]
                    myvalue = IncrementedValue(value[1], -mystep)
                    CopyDictCosmo[key] = [myvalue]
                    iterated = $(ReadMethod)(CopyDictCosmo, $(esc(Model)))
                    $(esc(VariedDict))["dvar_"*key*"_step_m_"*string(index)] = [iterated]
                end
            end
        end
        iterated = $(ReadMethod)($(esc(ReferenceDict)), $(esc(Model)))
        $(esc(VariedDict))["dvar_central_step_0"] = [iterated]
    end
end

"""
    CreateCosmologies(DictCosmo::Dict, steps::Array)

This function creates a dictionary containing all combinations required for forecasts.
"""
function CreateCosmologies(DictCosmo::Dict, CosmoModel::String, IADict::Dict, 
    IAModel::String, BiasDict::Dict, BiasModel::String, steps::Array)
    MyDictCosmo = Dict{String,Array{Any,1}}()
    @IterateParameters DictCosmo MyDictCosmo ReadCosmologyForecast steps CosmoModel

    MyDictIA = Dict{String,Array{Any,1}}()
    @IterateParameters IADict MyDictIA ReadIntrinsicAlignmentForecast steps IAModel

    MyDictBias = Dict{String,Array{Any,1}}()
    @IterateParameters BiasDict MyDictBias ReadBiasForecast steps BiasModel
    return MyDictCosmo, MyDictIA, MyDictBias
end

function CreateDirectoriesForecast!(Cosmologies::Dict, DictCosmo::Dict, IntrinsicAlignment::Dict, 
    DictIA::Dict, Bias::Dict, DictBias::Dict, path::String)
    mkdir(path)
    mkdir(path*"PowerSpectrum")
    CreateDirectoriesPmmCℓ(Cosmologies, path, "PowerSpectrum", true)
    mkdir(path*"Angular")
    CreateDirectoriesPmmCℓ(Cosmologies, path, "Angular", true)
    CreateDirectoriesPmmCℓ(IntrinsicAlignment, path, "Angular", false)
    CreateDirectoriesPmmCℓ(Bias, path, "Angular", false)
    mkdir(path*"Derivative/")
    CreateDirectoriesDerivative(DictCosmo, path)
    CreateDirectoriesDerivative(DictIA, path)
    CreateDirectoriesDerivative(DictBias, path)
end

function CreateDirectoriesPmmCℓ(DictionaryParams::Dict, PathFolder::String, NameFolder::String, 
    central::Bool)
    if central
        mkdir(PathFolder*"/"*NameFolder*"/"*"dvar_central_step_0")
    end
    for (key, value) in DictionaryParams
        if key != "dvar_central_step_0"
            mkdir(PathFolder*"/"*NameFolder*"/"*key)
        end
    end
end

function CreateDirectoriesDerivative(DictionaryParams::Dict, PathFolder::String)
    for (key, value) in DictionaryParams
        if value[2] == "present"
            mkdir(PathFolder*"Derivative/"*key)
        end
    end
end

function ForecastPowerSpectra!(Cosmologies::Dict, Path::String,
    CosmologicalGrid::CosmologicalGrid)
    for (key, value) in Cosmologies
        backgroundquantities = BackgroundQuantities(
        HZArray = zeros(length(CosmologicalGrid.ZArray)),
        rZArray=zeros(length(CosmologicalGrid.ZArray)))
        ComputeBackgroundQuantitiesGrid!(CosmologicalGrid,
        backgroundquantities, value[1])
        ClassyParams = Initializeclassy(value[1])
        powerspectrum = PowerSpectrum(PowerSpectrumLinArray =
        zeros(length(CosmologicalGrid.KArray), length(CosmologicalGrid.ZArray)),
        PowerSpectrumNonlinArray = zeros(length(CosmologicalGrid.KArray),
        length(CosmologicalGrid.ZArray)),
        InterpolatedPowerSpectrum = zeros(length(
        CosmologicalGrid.MultipolesArray), length(CosmologicalGrid.ZArray)))
        EvaluatePowerSpectrum!(ClassyParams, CosmologicalGrid, powerspectrum)
        WritePowerSpectrumBackground(powerspectrum, backgroundquantities,
        CosmologicalGrid, Path*key*"/p_mm")
        WriteCosmology!(value[1], Path*key)
    end
end

function InstantiateComputeWeightFunctionGrid(
    ConvolvedDensity::AbstractConvolvedDensity,
    w0waCDMCosmology::w0waCDMCosmology, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    GCWeightFunction::GCWeightFunction)
    ComputeBiasGrid!(CosmologicalGrid, GCWeightFunction, ConvolvedDensity)
    ComputeWeightFunctionGrid!(GCWeightFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
    return GCWeightFunction
end

function InstantiateComputeWeightFunctionGrid(
    ConvolvedDensity::AbstractConvolvedDensity,
    w0waCDMCosmology::w0waCDMCosmology,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    LensingFunction::WLWeightFunction, PathInput::String)
    ComputeLensingEfficiencyGrid!(LensingFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology, CustomLensingEfficiency())
    ComputeIntrinsicAlignmentGrid!(CosmologicalGrid, LensingFunction, ConvolvedDensity,
    BackgroundQuantities, w0waCDMCosmology, PathInput)
    ComputeWeightFunctionGrid!(LensingFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
    return LensingFunction
end

function Forecast∂Cℓ!(DictCosmo::Dict, IADict::Dict, BiasDict::Dict, PathInput::String,
    PathConfig::String, Steps::Array)
    ProbesDict = JSON.parsefile(PathConfig)
    CoefficientsArray = GetProbesArray(ProbesDict)
    for Coefficient in CoefficientsArray
        CentralCosmologyCL = ReadCℓ(
        PathInput*"/Angular/dvar_central_step_0/cl", Coefficient)
        CℓArray = zeros(size(CentralCosmologyCL.CℓArray, 1),
        size(CentralCosmologyCL.CℓArray, 2),
        size(CentralCosmologyCL.CℓArray, 3),
        2*length(Steps)+1)
        DerivativeArray = similar(CentralCosmologyCL.CℓArray)
        CℓArray[:, :, :, length(Steps) + 1] .=
        CentralCosmologyCL.CℓArray[:,:,:]
        stepvalues = zeros(2*length(Steps)+1)
        ∂cℓ = zeros(size(CℓArray, 1),
        size(CℓArray, 2), size(CℓArray, 3))
        #TODO: those three loops are quite similar. Macro?
        for (key, value) in DictCosmo
            if value[2] == "present"
                stepvalues[length(Steps) + 1] = value[1]
                for (index, mystep) in enumerate(Steps)
                    cℓminus = ReadCℓ(
                    PathInput*"/Angular/dvar_"*key*"_step_m_"*string(index)*"/cl",
                    Coefficient)
                    cℓplus = ReadCℓ(
                    PathInput*"/Angular/dvar_"*key*"_step_p_"*string(index)*"/cl",
                    Coefficient)
                    CℓArray[:, :, :, length(Steps) + 1 - index] .=
                    cℓminus.CℓArray
                    CℓArray[:, :, :, length(Steps) + 1 + index] .=
                    cℓplus.CℓArray
                    stepvalues[length(Steps) + 1 - index] =
                    IncrementedValue(value[1], -mystep)
                    stepvalues[length(Steps) + 1 + index] =
                    IncrementedValue(value[1],  mystep)
                end
                for idx_a in 1:size(CentralCosmologyCL.CℓArray, 2)
                    for idx_b in 1:size(CentralCosmologyCL.CℓArray, 3)
                        for idx_l in 1:size(CentralCosmologyCL.CℓArray, 1)
                            y = CℓArray[idx_l, idx_a, idx_b, :]
                            der = SteMDerivative(stepvalues, y)
                            ∂cℓ[idx_l, idx_a, idx_b] = der
                        end
                    end
                end
                Write∂Cℓ!(∂cℓ,
                PathInput*"/Derivative/"*key*"/"*key, Coefficient)
            end
        end
        for (key, value) in IADict
            if value[2] == "present"
                stepvalues[length(Steps) + 1] = value[1]
                for (index, mystep) in enumerate(Steps)
                    cℓminus = ReadCℓ(
                    PathInput*"/Angular/dvar_"*key*"_step_m_"*string(index)*"/cl",
                    Coefficient)
                    cℓplus = ReadCℓ(
                    PathInput*"/Angular/dvar_"*key*"_step_p_"*string(index)*"/cl",
                    Coefficient)
                    CℓArray[:, :, :, length(Steps) + 1 - index] .=
                    cℓminus.CℓArray
                    CℓArray[:, :, :, length(Steps) + 1 + index] .=
                    cℓplus.CℓArray
                    stepvalues[length(Steps) + 1 - index] =
                    IncrementedValue(value[1], -mystep)
                    stepvalues[length(Steps) + 1 + index] =
                    IncrementedValue(value[1],  mystep)
                end
                for idx_a in 1:size(CentralCosmologyCL.CℓArray, 2)
                    for idx_b in 1:size(CentralCosmologyCL.CℓArray, 3)
                        for idx_l in 1:size(CentralCosmologyCL.CℓArray, 1)
                            y = CℓArray[idx_l, idx_a, idx_b, :]
                            der = SteMDerivative(stepvalues, y)
                            ∂cℓ[idx_l, idx_a, idx_b] = der
                        end
                    end
                end
                Write∂Cℓ!(∂cℓ,
                PathInput*"/Derivative/"*key*"/"*key, Coefficient)
            end
        end
        for (key, value) in BiasDict
            if value[2] == "present"
                stepvalues[length(Steps) + 1] = value[1]
                for (index, mystep) in enumerate(Steps)
                    cℓminus = ReadCℓ(
                    PathInput*"/Angular/dvar_"*key*"_step_m_"*string(index)*"/cl",
                    Coefficient)
                    cℓplus = ReadCℓ(
                    PathInput*"/Angular/dvar_"*key*"_step_p_"*string(index)*"/cl",
                    Coefficient)
                    CℓArray[:, :, :, length(Steps) + 1 - index] .=
                    cℓminus.CℓArray
                    CℓArray[:, :, :, length(Steps) + 1 + index] .=
                    cℓplus.CℓArray
                    stepvalues[length(Steps) + 1 - index] =
                    IncrementedValue(value[1], -mystep)
                    stepvalues[length(Steps) + 1 + index] =
                    IncrementedValue(value[1],  mystep)
                end
                for idx_a in 1:size(CentralCosmologyCL.CℓArray, 2)
                    for idx_b in 1:size(CentralCosmologyCL.CℓArray, 3)
                        for idx_l in 1:size(CentralCosmologyCL.CℓArray, 1)
                            y = CℓArray[idx_l, idx_a, idx_b, :]
                            der = SteMDerivative(stepvalues, y)
                            ∂cℓ[idx_l, idx_a, idx_b] = der
                        end
                    end
                end
                Write∂Cℓ!(∂cℓ,
                PathInput*"/Derivative/"*key*"/"*key, Coefficient)
            end
        end
    end
end

#TODO: maybe we can create a macro?

function InstantiateWL(DictInput::Dict)
    LensingFunction = WLWeightFunction()
    return LensingFunction
end

function InstantiateWL(DictInput::Dict, IntrinsicAlignment::AbstractIntrinsicAlignment)
    LensingFunction = WLWeightFunction()
    LensingFunction.IntrinsicAlignmentModel = IntrinsicAlignment
    return LensingFunction
end

function InstantiateGC(DictInput::Dict)
    gcWeightFunction = GCWeightFunction()
    InstantiateBias(DictInput, gcWeightFunction)
    return gcWeightFunction
end

function InstantiateGC(DictInput::Dict, Bias::AbstractBias)
    gcWeightFunction = GCWeightFunction()
    gcWeightFunction.BiasKind = Bias
    return gcWeightFunction
end

function InstantiateBias(DictInput::Dict,
    GCWeightFunction::GCWeightFunction)
    if DictInput["PhotometricGalaxy"]["Bias"] == "PiecewiseBias"
        GCWeightFunction.BiasKind = PiecewiseBias()
    elseif DictInput["PhotometricGalaxy"]["Bias"] == "EuclidBias"
        GCWeightFunction.BiasKind = EuclidBias()
    else
        println("Bias must be correctly specified!")
    end
end

function InitializeProbes(DictInput::Dict,
    ConvolvedDensity::AbstractConvolvedDensity,
    w0waCDMCosmology::w0waCDMCosmology,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)
    DictProbes = Dict()
    if DictInput["Lensing"]["present"]
        LensingFunction = InstantiateWL(DictInput::Dict)
        LensingFunction = InstantiateComputeWeightFunctionGrid(ConvolvedDensity,
        w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities,
        LensingFunction, "../inputs/scaledmeanlum-E2Sa.txt")
        #TODO: remove the previous hardcoded line!
        push!(DictProbes, "Lensing" => LensingFunction)
    end
    if DictInput["PhotometricGalaxy"]["present"]
        GCWeightFunction = InstantiateGC(DictInput::Dict)
        GCWeightFunction =
        InstantiateComputeWeightFunctionGrid(ConvolvedDensity,
        w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities,
        GCWeightFunction)
        push!(DictProbes, "PhotometricGalaxy" => GCWeightFunction)
    end
    return DictProbes
end

function InitializeProbes(DictInput::Dict, ConvolvedDensity::AbstractConvolvedDensity,
    w0waCDMCosmology::w0waCDMCosmology, IntrinsicAlignment::AbstractIntrinsicAlignment,
    Bias::AbstractBias, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, PathInput::String)
    DictProbes = Dict()
    if DictInput["Lensing"]["present"]
        LensingFunction = InstantiateWL(DictInput::Dict)
        LensingFunction.IntrinsicAlignmentModel = IntrinsicAlignment
        LensingFunction = InstantiateComputeWeightFunctionGrid(ConvolvedDensity,
        w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities,
        LensingFunction, PathInput)
        push!(DictProbes, "Lensing" => LensingFunction)
    end
    if DictInput["PhotometricGalaxy"]["present"]
        GCWeightFunction = InstantiateGC(DictInput::Dict)
        GCWeightFunction.BiasKind = Bias
        GCWeightFunction = InstantiateComputeWeightFunctionGrid(ConvolvedDensity,
        w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities, GCWeightFunction)
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

function InitializeForecastCℓ(ProbesDict::Dict,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology,
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
                cℓ = Cℓ(CℓArray = 
                zeros(length(CosmologicalGrid.MultipolesArray),
                length(ProbesDict[key_A].WeightFunctionArray[:, 1]),
                length(ProbesDict[key_B].WeightFunctionArray[:, 1])))
                ComputeCℓ!(cℓ,
                ProbesDict[key_A], ProbesDict[key_B], BackgroundQuantities,
                w0waCDMCosmology, CosmologicalGrid, PowerSpectrum,
                CosmoCentral.CustomSimpson())
                WriteCℓ!(key_A*"_"*key_B,
                cℓ, PathOutput*key*"/cl")
                WriteCosmology!(w0waCDMCosmology, PathOutput*key)
            end
        end
    end
end

function InitializeForecastCℓ(ProbesDict::Dict,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology,
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
            #TODO is this if-else necessary?
            if key_B*"_"*key_A in CoefficientsArray
            else
                push!(CoefficientsArray, key_A*"_"*key_B)
                cℓ = Cℓ(CℓArray = zeros(length(CosmologicalGrid.MultipolesArray),
                length(ProbesDict[key_A].WeightFunctionArray[:, 1]),
                length(ProbesDict[key_B].WeightFunctionArray[:, 1])))
                ComputeCℓ!(cℓ, ProbesDict[key_A], ProbesDict[key_B], BackgroundQuantities,
                w0waCDMCosmology, CosmologicalGrid, PowerSpectrum,
                CosmoCentral.CustomSimpson())
                WriteCℓ!(key_A*"_"*key_B, cℓ, PathOutput*key*"/cl")
                WriteParameters!(CosmoDict, PathOutput*key)
            end
        end
    end
end

function ForecastCℓ!(Cosmologies::Dict, IntrinsicAlignment::Dict, Bias::Dict, 
    PathInputPmm::String, PathOutputCℓ::String, CosmologicalGrid::CosmologicalGrid,
    PathConfigCℓ::String)
    ProbesDict = JSON.parsefile(PathConfigCℓ)
    analyticaldensity = AnalitycalDensity()
    NormalizeAnalitycalDensity!(analyticaldensity)
    instrumentresponse = InstrumentResponse()
    convolveddensity = ConvolvedDensity(DensityGridArray =
    ones(10, length(CosmologicalGrid.ZArray)))
    NormalizeConvolvedDensity!(convolveddensity, analyticaldensity,
    instrumentresponse, CosmologicalGrid)
    ComputeConvolvedDensityGrid!(CosmologicalGrid, convolveddensity,
    analyticaldensity, instrumentresponse)
    #TODO: these three for loops are quite similar. The only difference is the cycled 
    #dictionary. Maybe we can create a macro?
    for (key, value) in Cosmologies
        w0waCDMCosmology = value[1]
        IA = IntrinsicAlignment["dvar_central_step_0"][1]
        bias = Bias["dvar_central_step_0"][1]
        PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
        ReadPowerSpectrumBackground(PathInputPmm*key*"/p_mm",
        CosmologicalGrid.MultipolesArray)
        CosmologicalGrid.MultipolesArray[1,1]
        DictProbes = InitializeProbes(ProbesDict, convolveddensity,
        w0waCDMCosmology, IA, bias, CosmologicalGrid, BackgroundQuantities, 
        "../inputs/scaledmeanlum-E2Sa.txt")
        ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
        InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
        BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
        InitializeForecastCℓ(DictProbes, BackgroundQuantities,
        w0waCDMCosmology, CosmologicalGrid, PowerSpectrum, PathOutputCℓ, key)
    end

    for (key, value) in IntrinsicAlignment
        if key != "dvar_central_step_0"
            w0waCDMCosmology = Cosmologies["dvar_central_step_0"][1]
            IA = value[1]
            bias = Bias["dvar_central_step_0"][1]
            PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
            ReadPowerSpectrumBackground(PathInputPmm*"dvar_central_step_0"*"/p_mm",
            CosmologicalGrid.MultipolesArray)
            CosmologicalGrid.MultipolesArray[1,1]
            DictProbes = InitializeProbes(ProbesDict, convolveddensity,
            w0waCDMCosmology, IA, bias, CosmologicalGrid, BackgroundQuantities, 
            "../inputs/scaledmeanlum-E2Sa.txt")
            ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
            InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
            BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
            InitializeForecastCℓ(DictProbes, BackgroundQuantities,
            w0waCDMCosmology, CosmologicalGrid, PowerSpectrum, PathOutputCℓ, key)
        end
    end

    for (key, value) in Bias
        if key != "dvar_central_step_0"
            w0waCDMCosmology = Cosmologies["dvar_central_step_0"][1]
            IA = IntrinsicAlignment["dvar_central_step_0"][1]
            bias = value[1]
            PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
            ReadPowerSpectrumBackground(PathInputPmm*"dvar_central_step_0"*"/p_mm",
            CosmologicalGrid.MultipolesArray)
            CosmologicalGrid.MultipolesArray[1,1]
            DictProbes = InitializeProbes(ProbesDict, convolveddensity,
            w0waCDMCosmology, IA, bias, CosmologicalGrid, BackgroundQuantities, 
            "../inputs/scaledmeanlum-E2Sa.txt")
            ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
            InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
            BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
            InitializeForecastCℓ(DictProbes, BackgroundQuantities,
            w0waCDMCosmology, CosmologicalGrid, PowerSpectrum, PathOutputCℓ, key)
        end
    end
end
