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
        χZArray=zeros(length(CosmologicalGrid.ZArray)))
        ComputeBackgroundQuantitiesGrid!(CosmologicalGrid,
        backgroundquantities, value[1])
        ClassyParams = Initializeclassy(value[1])
        powerspectrum = PowerSpectrum(PowerSpectrumLinArray =
        zeros(length(CosmologicalGrid.KArray), length(CosmologicalGrid.ZArray)),
        PowerSpectrumNonlinArray = zeros(length(CosmologicalGrid.KArray),
        length(CosmologicalGrid.ZArray)),
        InterpolatedPowerSpectrum = zeros(length(
        CosmologicalGrid.ℓBinCenters), length(CosmologicalGrid.ZArray)))
        EvaluatePowerSpectrum!(ClassyParams, CosmologicalGrid, powerspectrum)
        WritePowerSpectrumBackground(powerspectrum, backgroundquantities,
        CosmologicalGrid, Path*key*"/p_mm")
        WriteCosmology!(value[1], Path*key)
    end
end

function InstantiateComputeWeightFunctionGrid(
    ConvolvedDensity::AbstractConvolvedDensity,
    Cosmology::AbstractCosmology, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    GCWeightFunction::GCWeightFunction)
    ComputeBiasGrid!(CosmologicalGrid, GCWeightFunction, ConvolvedDensity)
    ComputeWeightFunctionGrid!(GCWeightFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, Cosmology)
    return GCWeightFunction
end

function InstantiateComputeWeightFunctionGrid(
    ConvolvedDensity::AbstractConvolvedDensity,
    Cosmology::AbstractCosmology,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    LensingFunction::WLWeightFunction)
    ComputeLensingEfficiencyGrid!(LensingFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, Cosmology, CustomLensingEfficiency())
    ComputeIntrinsicAlignmentGrid!(CosmologicalGrid, LensingFunction, ConvolvedDensity,
    BackgroundQuantities, Cosmology)
    ComputeWeightFunctionGrid!(LensingFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, Cosmology)
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
    Cosmology::AbstractCosmology,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities)
    DictProbes = Dict()
    if DictInput["Lensing"]["present"]
        LensingFunction = InstantiateWL(DictInput::Dict)
        LensingFunction = InstantiateComputeWeightFunctionGrid(ConvolvedDensity,
        Cosmology, CosmologicalGrid, BackgroundQuantities,
        LensingFunction)
        push!(DictProbes, "Lensing" => LensingFunction)
    end
    if DictInput["PhotometricGalaxy"]["present"]
        GCWeightFunction = InstantiateGC(DictInput::Dict)
        GCWeightFunction =
        InstantiateComputeWeightFunctionGrid(ConvolvedDensity,
        Cosmology, CosmologicalGrid, BackgroundQuantities,
        GCWeightFunction)
        push!(DictProbes, "PhotometricGalaxy" => GCWeightFunction)
    end
    return DictProbes
end

function InitializeProbes(DictInput::Dict, ConvolvedDensity::AbstractConvolvedDensity,
    Cosmology::AbstractCosmology, IntrinsicAlignment::AbstractIntrinsicAlignment,
    Bias::AbstractBias, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, PathInput::String)
    DictProbes = Dict()
    if DictInput["Lensing"]["present"]
        LensingFunction = InstantiateWL(DictInput::Dict)
        LensingFunction.IntrinsicAlignmentModel = IntrinsicAlignment
        LensingFunction = InstantiateComputeWeightFunctionGrid(ConvolvedDensity,
        Cosmology, CosmologicalGrid, BackgroundQuantities,
        LensingFunction)
        push!(DictProbes, "Lensing" => LensingFunction)
    end
    if DictInput["PhotometricGalaxy"]["present"]
        GCWeightFunction = InstantiateGC(DictInput::Dict)
        GCWeightFunction.BiasKind = Bias
        GCWeightFunction = InstantiateComputeWeightFunctionGrid(ConvolvedDensity,
        Cosmology, CosmologicalGrid, BackgroundQuantities, GCWeightFunction)
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
    Cosmology::AbstractCosmology,
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
                zeros(length(CosmologicalGrid.ℓBinCenters),
                length(ProbesDict[key_A].WeightFunctionArray[:, 1]),
                length(ProbesDict[key_B].WeightFunctionArray[:, 1])))
                ComputeCℓ!(cℓ,
                ProbesDict[key_A], ProbesDict[key_B], BackgroundQuantities,
                Cosmology, CosmologicalGrid, PowerSpectrum,
                CosmoCentral.CustomSimpson())
                WriteCℓ!(key_A*"_"*key_B,
                cℓ, PathOutput*key*"/cl")
                WriteCosmology!(Cosmology, PathOutput*key)
            end
        end
        WriteWeightFunctions!(key_A, ProbesDict[key_A], PathOutput*key*"/cl")
    end
    WriteCosmologicalGrid!(PathOutput*key*"/cl", CosmologicalGrid)
end

function InitializeForecastCℓ(ProbesDict::Dict,
    BackgroundQuantities::BackgroundQuantities,
    Cosmology::AbstractCosmology,
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
                cℓ = Cℓ(CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters),
                length(ProbesDict[key_A].WeightFunctionArray[:, 1]),
                length(ProbesDict[key_B].WeightFunctionArray[:, 1])))
                ComputeCℓ!(cℓ, ProbesDict[key_A], ProbesDict[key_B], BackgroundQuantities,
                Cosmology, CosmologicalGrid, PowerSpectrum,
                CosmoCentral.CustomSimpson())
                WriteCℓ!(key_A*"_"*key_B, cℓ, PathOutput*key*"/cl")
                WriteParameters!(CosmoDict, PathOutput*key)
            end
        end
    end
end

function ForecastCℓ!(Cosmologies::Dict, IntrinsicAlignment::Dict, Bias::Dict, 
    PathInputPmm::String, PathOutputCℓ::String, CosmologicalGrid::CosmologicalGrid,
    PathConfigCℓ::String, PathInputCℓ::String)
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
    #dictionary. Maybe we can create a macro? The only difficulty is that we are iterating
    #over different dictionaries. maybe we can write some if in the dficitonaries to decide
    #which expr return?
    for (key, value) in Cosmologies
        Cosmology = value[1]
        IA = IntrinsicAlignment["dvar_central_step_0"][1]
        bias = Bias["dvar_central_step_0"][1]
        PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
        ReadPowerSpectrumBackground(PathInputPmm*key*"/p_mm",
        CosmologicalGrid.ℓBinCenters, CosmologicalGrid.ℓBinWidths)
        ExtractGrowthFactor!(BackgroundQuantities, PowerSpectrum)
        CosmologicalGrid.ℓBinCenters[1,1] #TODO wtf represents this???
        DictProbes = InitializeProbes(ProbesDict, convolveddensity, Cosmology, IA,
        bias, CosmologicalGrid, BackgroundQuantities, PathInputCℓ)
        ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
        InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid, BackgroundQuantities,
        PowerSpectrum, CosmoCentral.BSplineCubic())
        InitializeForecastCℓ(DictProbes, BackgroundQuantities, Cosmology,
        CosmologicalGrid, PowerSpectrum, PathOutputCℓ, key)
    end

    for (key, value) in IntrinsicAlignment
        if key != "dvar_central_step_0"
            Cosmology = Cosmologies["dvar_central_step_0"][1]
            IA = value[1]
            bias = Bias["dvar_central_step_0"][1]
            PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
            ReadPowerSpectrumBackground(PathInputPmm*"dvar_central_step_0"*"/p_mm",
            CosmologicalGrid.ℓBinCenters, CosmologicalGrid.ℓBinWidths)
            ExtractGrowthFactor!(BackgroundQuantities, PowerSpectrum)
            CosmologicalGrid.ℓBinCenters[1,1]
            DictProbes = InitializeProbes(ProbesDict, convolveddensity,
            Cosmology, IA, bias, CosmologicalGrid, BackgroundQuantities, 
            PathInputCℓ)
            ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
            InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
            BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
            InitializeForecastCℓ(DictProbes, BackgroundQuantities,
            Cosmology, CosmologicalGrid, PowerSpectrum, PathOutputCℓ, key)
        end
    end

    for (key, value) in Bias
        if key != "dvar_central_step_0"
            Cosmology = Cosmologies["dvar_central_step_0"][1]
            IA = IntrinsicAlignment["dvar_central_step_0"][1]
            bias = value[1]
            PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
            ReadPowerSpectrumBackground(PathInputPmm*"dvar_central_step_0"*"/p_mm",
            CosmologicalGrid.ℓBinCenters, CosmologicalGrid.ℓBinWidths)
            ExtractGrowthFactor!(BackgroundQuantities, PowerSpectrum)
            CosmologicalGrid.ℓBinCenters[1,1]
            DictProbes = InitializeProbes(ProbesDict, convolveddensity,
            Cosmology, IA, bias, CosmologicalGrid, BackgroundQuantities, 
            PathInputCℓ)
            ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
            InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
            BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
            InitializeForecastCℓ(DictProbes, BackgroundQuantities,
            Cosmology, CosmologicalGrid, PowerSpectrum, PathOutputCℓ, key)
        end
    end
end

function ExtractVariedParameters(InputList::Vector{Dict{String, Vector{Any}}})
    OutputList = []
    for myDict in InputList
        for (key, value) in myDict
            if value[2] == "present"
                push!(OutputList, key)
            end
        end
    end
    return OutputList
end

"""
    ForecastFisherαβ(PathCentralCℓ::String, Path∂Cℓ::String,
    InputList::Vector{Dict{String, Vector{Any}}}, CosmologicalGrid::CosmologicalGrid)

This function evaluate the Fisher Matrix according to the following formula:
```math
F_{\\alpha \\beta}=\\sum_{\\ell=\\ell_{\\min }}^{\\ell_{\\max }} \\operatorname{Tr}\\left
\\{[\\boldsymbol{\\Sigma}(\\ell)]^{-1} \\frac{\\partial \\mathbf{C}(\\ell)}{\\partial 
\\alpha}[\\boldsymbol{\\Sigma}(\\ell)]^{-1} \\frac{\\partial \\mathbf{C}(\\ell)}{\\partial
\\beta}\\right\\},
```
where ``\\alpha`` and ``\\beta`` are parameters of the Fisher Matrix, ``\\ell`` are the
multipoles, ``\\Sigma`` is the Covariance Matrix of the  Field Approach and
``\\frac{\\partial \\mathbf{C}(\\ell)}{\\partial \\alpha}`` is the Matrix of the derivatives
of the ``C_\\ell`` wrt parameter ``\\alpha``.
"""
function ForecastFisherαβ(PathCentralCℓ::String, Path∂Cℓ::String,
    InputList::Vector{Dict{String, Vector{Any}}}, CosmologicalGrid::CosmologicalGrid)
    Fisher = Fisherαβ()
    VariedParameters = ExtractVariedParameters(InputList)
    Fisher.FisherMatrix = zeros(length(VariedParameters), length(VariedParameters))
    Fisher.ParametersList = VariedParameters
    Fisher.SelectedParametersList = VariedParameters
    AnalitycalDensity = CosmoCentral.AnalitycalDensity()
    NormalizeAnalitycalDensity!(AnalitycalDensity)
    InstrumentResponse = CosmoCentral.InstrumentResponse()
    ConvolvedDensity = CosmoCentral.ConvolvedDensity(DensityGridArray =
    ones(10, length(CosmologicalGrid.ZArray)))
    NormalizeConvolvedDensity!(ConvolvedDensity, AnalitycalDensity, InstrumentResponse,
    CosmologicalGrid)
    ComputeConvolvedDensityGrid!(CosmologicalGrid, ConvolvedDensity, AnalitycalDensity,
    InstrumentResponse)
    ComputeSurfaceDensityBins!(ConvolvedDensity, AnalitycalDensity)

    Cℓ = ReadCℓ(PathCentralCℓ, "Lensing_Lensing")
    #TODO now only LL, but this need to be more flexible...maybe list with probes?
    Cov = InstantiateEvaluateCovariance(Cℓ, ConvolvedDensity, CosmologicalGrid, "Lensing",
    "Lensing")
    EvaluateFisherMatrix!(VariedParameters, Fisher, Path∂Cℓ, Cov)
    SelectMatrixAndMarginalize!(VariedParameters, Fisher)
    return Fisher
end


"""
    ForecastFisherαβ(PathCentralCℓ::String, Path∂Cℓ::String,
    InputList::Vector{Dict{String, Vector{Any}}}, CosmologicalGrid::CosmologicalGrid,
    ciccio::String)

This function evaluate the Fisher Matrix according to the following formula:
```math
F_{\\alpha \\beta}=\\sum_{\\ell=\\ell_{\\min }}^{\\ell_{\\max }} \\operatorname{vecp}
\\left(\\frac{\\partial \\mathbf{C}(\\ell)}{\\partial \\alpha}\\right)^{T}
\\left(\\boldsymbol{\\Xi}(\\ell)\\right)^{-1} \\operatorname{vecp}\\left(\\frac{\\partial
\\mathbf{C}(\\ell)} {\\partial \\beta}\\right),
```
where ``\\alpha`` and ``\\beta`` are parameters of the Fisher Matrix, ``\\ell`` are the
multipoles, ``\\Xi`` is the Covariance Matrix of the  Field Approach
[`CℓCovariance`](@ref) and ``\\frac{\\partial \\mathbf{C}(\\ell)}{\\partial \\alpha}``
is the Matrix of the derivatives of the ``C_\\ell`` wrt parameter ``\\alpha``.
"""
function ForecastFisherαβ(PathCentralCℓ::String, Path∂Cℓ::String,
    InputList::Vector{Dict{String, Vector{Any}}}, CosmologicalGrid::CosmologicalGrid,
    ciccio::String)
    Fisher = Fisherαβ()
    VariedParameters = ExtractVariedParameters(InputList)
    Fisher.FisherMatrix = zeros(length(VariedParameters), length(VariedParameters))
    Fisher.ParametersList = VariedParameters
    Fisher.SelectedParametersList = VariedParameters
    #here we instantiate the density again to evaluate the noise. Maybe it could be better
    #if we evaluated the density once for all, we passed it to Cℓ evaluator and then to
    #Fisher, in order to be more safe.
    AnalitycalDensity = CosmoCentral.AnalitycalDensity()
    NormalizeAnalitycalDensity!(AnalitycalDensity)
    InstrumentResponse = CosmoCentral.InstrumentResponse()
    ConvolvedDensity = CosmoCentral.ConvolvedDensity(DensityGridArray =
    ones(10, length(CosmologicalGrid.ZArray)))
    NormalizeConvolvedDensity!(ConvolvedDensity, AnalitycalDensity, InstrumentResponse,
    CosmologicalGrid)
    ComputeConvolvedDensityGrid!(CosmologicalGrid, ConvolvedDensity, AnalitycalDensity,
    InstrumentResponse)
    ComputeSurfaceDensityBins!(ConvolvedDensity, AnalitycalDensity)

    Cℓ = ReadCℓ(PathCentralCℓ, "Lensing_Lensing")
    #TODO now only LL, but this need to be more flexible...maybe list with probes?
    Covaₗₘ = InstantiateEvaluateCovariance(Cℓ, ConvolvedDensity, CosmologicalGrid, "Lensing",
    "Lensing")
    CovCℓ = InstantiateEvaluateCovariance(Covaₗₘ)
    EvaluateFisherMatrix!(VariedParameters, Fisher, Path∂Cℓ, CovCℓ)
    SelectMatrixAndMarginalize!(VariedParameters, Fisher)
    return Fisher
end

function EvaluateFisherMatrix!(VariedParameters::Vector{}, Fisher::AbstractFisher,
    Path∂Cℓ::String, Cov::AbstractCovariance)
    for (indexα, Parα) in enumerate(VariedParameters)
        ∂Cℓα = Read∂Cℓ(Path∂Cℓ*"/"*Parα*"/"*Parα, "Lensing_Lensing")
        for (indexβ, Parβ) in enumerate(VariedParameters)
            ∂Cℓβ = Read∂Cℓ(Path∂Cℓ*"/"*Parβ*"/"*Parβ, "Lensing_Lensing")
            EvaluateFisherMatrixElement!(Fisher, Cov, ∂Cℓα, ∂Cℓβ, Parα, Parβ)
            Fisher.FisherMatrix[indexα, indexβ] = Fisher.FisherDict[Parα*"_"*Parβ]
        end
    end
end

function SelectMatrixAndMarginalize!(VariedParameters::Vector{}, Fisher::AbstractFisher)
    for (idx, Par) in enumerate(reverse(VariedParameters))
        reverse_index = length(VariedParameters)-idx+1
        if Fisher.FisherMatrix[reverse_index, reverse_index] == 0
            Fisher.FisherMatrix = Fisher.FisherMatrix[1:end .!= reverse_index,
            1:end .!= reverse_index]
            Fisher.SelectedParametersList =
            Fisher.SelectedParametersList[1:end .!= reverse_index]
        end
    end
    Fisher.CorrelationMatrix = inv(Fisher.FisherMatrix)
    for (idxα, Parα) in enumerate(Fisher.SelectedParametersList)
        Fisher.MarginalizedErrors[Parα] = sqrt(Fisher.CorrelationMatrix[idxα, idxα])
    end
end

#Here we try to implement the NEW forecaster, more manageable and flexible
#Before removing the old one, we will check that everything works

function ReadPresentDict(dict::Dict)
    newdict = Dict()
    for (key,value) in dict
        if key == "model"
            newdict["model"] = value
        else
            newdict[key] = value["value"]
        end
    end
    return newdict
end

function Iterator!(ContainerDict::Dict, BaseDict::Dict, SteMSteps::Array)
    for (key,value) in BaseDict
        if key != "model"
            if value["free"]
                for (index, mystep) in enumerate(SteMSteps)
                    ContainerDict["dvar_"*key*"_step_p_"*string(index)] =
                    deepcopy(ReadPresentDict(BaseDict))
                    ContainerDict["dvar_"*key*"_step_p_"*string(index)][key] = 
                    IncrementedValue(BaseDict[key]["value"], mystep)
                    ContainerDict["dvar_"*key*"_step_m_"*string(index)] =
                    deepcopy(ReadPresentDict(BaseDict))
                    ContainerDict["dvar_"*key*"_step_m_"*string(index)][key] =
                    IncrementedValue(BaseDict[key]["value"], -mystep)
                end
            end
        end
    end
end

function ReadInputForecast(InputDict::Dict, cosmogrid, SteMSteps)
    ProbesDict = Dict()
    for (key,value) in InputDict
        ProbesDict[key] = IterateProbe(value, cosmogrid, SteMSteps)
    end
    return ProbesDict
end

function CreateCosmology(CosmoDict::Dict)
    if CosmoDict["model"] == "Flatw0waCDM"
        Cosmo = CreateFlatw0waCDM((CosmoDict))
    else
        error("No Cosmology!")
    end
    return Cosmo
end

function IterateCosmologies(CosmoDict::Dict, SteMSteps::Array)
    CosmoVariedDict = Dict()
    CosmoVariedDict["dvar_central_step_0"] = (ReadPresentDict(CosmoDict))
    Iterator!(CosmoVariedDict, CosmoDict, SteMSteps)
    return CosmoVariedDict
end

function CreateFlatw0waCDM(CosmoDict::Dict)
    Cosmo = CosmoCentral.Flatw0waCDMCosmology()
    Cosmo.w0 = CosmoDict["w0"]
    Cosmo.wa = CosmoDict["wa"]
    Cosmo.Mν = CosmoDict["Mν"]
    Cosmo.H0 = CosmoDict["H0"]
    Cosmo.ΩM = CosmoDict["ΩM"]
    Cosmo.ΩB = CosmoDict["ΩB"]
    Cosmo.ns = CosmoDict["ns"]
    Cosmo.σ8 = CosmoDict["σ8"]
    return Cosmo
end

function ReadInputProbesForecast(InputDict::Dict, cosmogrid::CosmologicalGrid,
    SteMSteps::Array)
    ProbesDict = Dict()
    for (key,value) in InputDict
        ProbesDict[key] = IterateProbe(value, cosmogrid, SteMSteps)
    end
    return ProbesDict
end

abstract type AbstractContainer end

@kwdef mutable struct ForecastContainer <: AbstractContainer
    CosmoDict::Dict = Dict()
    ProbesDict::Dict = Dict()
    VariedParsList::Array = []
end

function InitializeForecastContainer(CosmoDict::Dict, ProbesDict::Dict,
    cosmogrid::CosmologicalGrid, SteMSteps::Array)
    forcontainer = ForecastContainer()
    forcontainer.CosmoDict = IterateCosmologies(CosmoDict, SteMSteps)
    forcontainer.ProbesDict = ReadInputProbesForecast(ProbesDict, cosmogrid, SteMSteps)
    SelectVariedParameters!(forcontainer, CosmoDict, ProbesDict)
    return forcontainer
end

function SelectSubDictProbe(ProbeDict::Dict)
    if ProbeDict["probe"] == "WLProbe"
        return ProbeDict["IntrinsicAlignment"]
    elseif ProbeDict["probe"] == "GCProbe"
        return ProbeDict["Bias"]
    else
        error("Not selected probe")
    end
end

function SelectVariedParameters!(forcontainer::ForecastContainer, CosmoDict::Dict, ProbesDict::Dict)
    varied_parslist = []
    for (key, value) in CosmoDict
        if key != "model"
            if value["free"] == true
                push!(varied_parslist, (key, value["value"]))
            end
        end
    end
    for (keyprobe, valueprobe) in ProbesDict
        subdict = SelectSubDictProbe(valueprobe)
        for (key, value) in subdict
            if key != "model"
                if value["free"] == true
                    push!(varied_parslist, (key, value["value"]))
                end
            end
        end
    end
    forcontainer.VariedParsList = varied_parslist
end

function CreateDirectoriesForecast!(forcontainer::ForecastContainer, path::String)
    mkdir(path)
    mkdir(path*"PowerSpectrum")
    mkdir(path*"Angular")
    mkdir(path*"Derivative")
    for (key,value) in forcontainer.CosmoDict
        mkdir(path*"PowerSpectrum/"*key)
        mkdir(path*"Angular/"*key)
    end
    for (keyprobe,valueprobe) in forcontainer.ProbesDict
        for (keyvar,valuevar) in valueprobe
            if keyvar != "dvar_central_step_0"
                mkdir(path*"Angular/"*keyvar)
            end
        end
    end
    for (key, value) in forcontainer.VariedParsList
        mkdir(path*"Derivative/"*key)
    end
end

function ForecastPowerSpectra!(forcontainer::ForecastContainer, Path::String,
    CosmologicalGrid::CosmologicalGrid)
    for (key, value) in forcontainer.CosmoDict
        cosmology = CreateCosmology(value)
        backgroundquantities = BackgroundQuantities(
        HZArray = zeros(length(CosmologicalGrid.ZArray)),
        χZArray=zeros(length(CosmologicalGrid.ZArray)))
        ComputeBackgroundQuantitiesGrid!(CosmologicalGrid,
        backgroundquantities, cosmology)
        ClassyParams = Initializeclassy(cosmology)
        powerspectrum = PowerSpectrum(PowerSpectrumLinArray =
        zeros(length(CosmologicalGrid.KArray), length(CosmologicalGrid.ZArray)),
        PowerSpectrumNonlinArray = zeros(length(CosmologicalGrid.KArray),
        length(CosmologicalGrid.ZArray)),
        InterpolatedPowerSpectrum = zeros(length(
        CosmologicalGrid.ℓBinCenters), length(CosmologicalGrid.ZArray)))
        EvaluatePowerSpectrum!(ClassyParams, CosmologicalGrid, powerspectrum)
        WritePowerSpectrumBackground(powerspectrum, backgroundquantities,
        CosmologicalGrid, Path*key*"/p_mm")
        WriteCosmology!(cosmology, Path*key)
    end
end

function CreateProbesAndObservablesArray(forcontainer::ForecastContainer)
    probes_array = []
    for (keyprobe,valueprobe) in forcontainer.ProbesDict
        push!(probes_array, keyprobe)
    end
    sort!(probes_array)
    observables_array = []
    for probea in probes_array
        for probeb in probes_array
            if probeb*"_"*probea ∉ observables_array
                push!(observables_array, probea*"_"*probeb)
            end
        end
    end
    return probes_array, observables_array
end

function ForecastCℓ!(forcontainer::ForecastContainer, cosmogrid::CosmologicalGrid,
    PathInputPmm::String, PathOutputCℓ::String)
    probes_array, observable_array = CreateProbesAndObservablesArray(forcontainer)
    weightdict = Dict()
    for (keycosmo,valuecosmo) in forcontainer.CosmoDict
        cosmology = ReadCosmology(valuecosmo)
        PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
        CosmoCentral.ReadPowerSpectrumBackground(PathInputPmm*keycosmo*"/p_mm",
        cosmogrid.ℓBinCenters, cosmogrid.ℓBinWidths)
        CosmoCentral.ExtractGrowthFactor!(BackgroundQuantities, PowerSpectrum)
        ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
        InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid, BackgroundQuantities,
        PowerSpectrum, BSplineCubic())
        for (keyprobe, valueprobe) in forcontainer.ProbesDict
            for probe in probes_array
                weightdict[probe] = CreateAndEvaluateWeightFunction(
                forcontainer.ProbesDict[probe]["dvar_central_step_0"], CosmologicalGrid,
                BackgroundQuantities, cosmology)
            end
        end
        used_probes = []
        for probea in probes_array
            WriteWeightFunctions!(probea, weightdict[probea],
                PathOutputCℓ*keycosmo*"/cl")
            for probeb in probes_array
                if probeb*"_"*probea ∉ used_probes
                    push!(used_probes, probea*"_"*probeb)
                    cℓ = Cℓ(CℓArray = 
                    zeros(length(CosmologicalGrid.ℓBinCenters),
                    length(weightdict[probea].WeightFunctionArray[:, 1]),
                    length(weightdict[probeb].WeightFunctionArray[:, 1])))
                    ComputeCℓ!(cℓ, weightdict[probea], weightdict[probeb],
                    BackgroundQuantities, cosmology, CosmologicalGrid, PowerSpectrum,
                    CustomSimpson())
                    WriteCℓ!(probea*"_"*probeb, cℓ, PathOutputCℓ*keycosmo*"/cl")
                    WriteCosmology!(cosmology, PathOutputCℓ*keycosmo)
                end
            end
        end
    end
    keycosmo = "dvar_central_step_0"
    cosmology = ReadCosmology(forcontainer.CosmoDict[keycosmo])
    PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
    CosmoCentral.ReadPowerSpectrumBackground(PathInputPmm*keycosmo*"/p_mm",
    cosmogrid.ℓBinCenters, cosmogrid.ℓBinWidths)
    CosmoCentral.ExtractGrowthFactor!(BackgroundQuantities, PowerSpectrum)
    ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
    InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid, BackgroundQuantities,
    PowerSpectrum, BSplineCubic())
    for probea in probes_array
        for (keyvarieda,valuevarieda) in forcontainer.ProbesDict[probea]
            if keyvarieda != "dvar_central_step_0"
                weightdict[probea] = CreateAndEvaluateWeightFunction(valuevarieda,
                CosmologicalGrid, BackgroundQuantities, cosmology)
                WriteWeightFunctions!(probea, weightdict[probea],
                PathOutputCℓ*keyvarieda*"/cl")
                for probeb in probes_array
                    if probeb== probea
                        cℓ = Cℓ(CℓArray = 
                        zeros(length(CosmologicalGrid.ℓBinCenters),
                        length(weightdict[probea].WeightFunctionArray[:, 1]),
                        length(weightdict[probea].WeightFunctionArray[:, 1])))
                        ComputeCℓ!(cℓ, weightdict[probea], weightdict[probea],
                        BackgroundQuantities, cosmology, CosmologicalGrid, PowerSpectrum,
                        CustomSimpson())
                        WriteCℓ!(probea*"_"*probea, cℓ, PathOutputCℓ*keyvarieda*"/cl")
                        WriteCosmology!(cosmology, PathOutputCℓ*keyvarieda)

                    else
                        weightdict[probeb] = CreateAndEvaluateWeightFunction(
                        forcontainer.ProbesDict[probeb]["dvar_central_step_0"],
                        CosmologicalGrid, BackgroundQuantities, cosmology)
                        WriteWeightFunctions!(probeb, weightdict[probeb],
                        PathOutputCℓ*keyvarieda*"/cl")
                        cℓ = Cℓ(CℓArray = 
                        zeros(length(CosmologicalGrid.ℓBinCenters),
                        length(weightdict[probea].WeightFunctionArray[:, 1]),
                        length(weightdict[probeb].WeightFunctionArray[:, 1])))
                        ComputeCℓ!(cℓ, weightdict[probea], weightdict[probeb],
                        BackgroundQuantities, cosmology, CosmologicalGrid, PowerSpectrum,
                        CustomSimpson())
                        if probea < probeb
                            WriteCℓ!(probea*"_"*probeb, cℓ, PathOutputCℓ*keyvarieda*"/cl")
                        else
                            WriteCℓ!(probeb*"_"*probea, cℓ, PathOutputCℓ*keyvarieda*"/cl")
                        end
                        WriteCosmology!(cosmology, PathOutputCℓ*keyvarieda)
                        cℓ = Cℓ(CℓArray = 
                        zeros(length(CosmologicalGrid.ℓBinCenters),
                        length(weightdict[probea].WeightFunctionArray[:, 1]),
                        length(weightdict[probeb].WeightFunctionArray[:, 1])))
                        ComputeCℓ!(cℓ, weightdict[probeb], weightdict[probeb],
                        BackgroundQuantities, cosmology, CosmologicalGrid, PowerSpectrum,
                        CustomSimpson())
                        WriteCℓ!(probeb*"_"*probeb, cℓ, PathOutputCℓ*keyvarieda*"/cl")
                    end
                end
            end
        end
    end
end

function Forecast∂Cℓ!(forcontainer::ForecastContainer, PathInput::String,
    PathConfig::String, Steps::Array)
    #ProbesDict = JSON.parsefile(PathConfig)
    #CoefficientsArray = GetProbesArray(ProbesDict)
    probes_array, observable_array = CreateProbesAndObservablesArray(forcontainer)
    for observable in observable_array
        CentralCosmologyCL = ReadCℓ(
        PathInput*"/Angular/dvar_central_step_0/cl", observable)
        CℓArray = zeros(size(CentralCosmologyCL.CℓArray, 1),
        size(CentralCosmologyCL.CℓArray, 2),
        size(CentralCosmologyCL.CℓArray, 3),
        2*length(Steps)+1)
        DerivativeArray = zeros(size(CentralCosmologyCL.CℓArray))
        CℓArray[:, :, :, length(Steps) + 1] .=
        CentralCosmologyCL.CℓArray[:,:,:]
        stepvalues = zeros(2*length(Steps)+1)
        ∂cℓ = zeros(size(CℓArray, 1),
        size(CℓArray, 2), size(CℓArray, 3))
        for (key, value) in forcontainer.VariedParsList
            stepvalues[length(Steps) + 1] = value
            for (index, mystep) in enumerate(Steps)
                cℓminus = ReadCℓ(
                PathInput*"/Angular/dvar_"*key*"_step_m_"*string(index)*"/cl",
                observable)
                cℓplus = ReadCℓ(
                PathInput*"/Angular/dvar_"*key*"_step_p_"*string(index)*"/cl",
                observable)
                CℓArray[:, :, :, length(Steps) + 1 - index] .=
                cℓminus.CℓArray
                CℓArray[:, :, :, length(Steps) + 1 + index] .=
                cℓplus.CℓArray
                stepvalues[length(Steps) + 1 - index] =
                IncrementedValue(value, -mystep)
                stepvalues[length(Steps) + 1 + index] =
                IncrementedValue(value,  mystep)
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
            PathInput*"/Derivative/"*key*"/"*key, observable)
        end
    end
end

function ExtractVariedParameters(forcontainer::ForecastContainer)
    OutputList = []
    for (key, value) in forcontainer.VariedParsList
        push!(OutputList, key)
    end
    return OutputList
end

function ForecastFisherαβ(forcontainer::ForecastContainer ,PathCentralCℓ::String, Path∂Cℓ::String,
    CosmologicalGrid::CosmologicalGrid, ciccio::String)
    Fisher = Fisherαβ()
    VariedParameters = ExtractVariedParameters(forcontainer)
    Fisher.FisherMatrix = zeros(length(VariedParameters), length(VariedParameters))
    Fisher.ParametersList = VariedParameters
    Fisher.SelectedParametersList = VariedParameters
    #here we instantiate the density again to evaluate the noise. Maybe it could be better
    #if we evaluated the density once for all, we passed it to Cℓ evaluator and then to
    #Fisher, in order to be more safe.
    AnalitycalDensity = CosmoCentral.AnalitycalDensity()
    NormalizeAnalitycalDensity!(AnalitycalDensity)
    InstrumentResponse = CosmoCentral.InstrumentResponse()
    ConvolvedDensity = CosmoCentral.ConvolvedDensity(DensityGridArray =
    ones(10, length(CosmologicalGrid.ZArray)))
    NormalizeConvolvedDensity!(ConvolvedDensity, AnalitycalDensity, InstrumentResponse,
    CosmologicalGrid)
    ComputeConvolvedDensityGrid!(CosmologicalGrid, ConvolvedDensity, AnalitycalDensity,
    InstrumentResponse)
    ComputeSurfaceDensityBins!(ConvolvedDensity, AnalitycalDensity)

    Cℓ = ReadCℓ(PathCentralCℓ, "Lensing_Lensing")
    #TODO now only LL, but this need to be more flexible...maybe list with probes?
    Covaₗₘ = InstantiateEvaluateCovariance(Cℓ, ConvolvedDensity, CosmologicalGrid, "Lensing",
    "Lensing")
    CovCℓ = InstantiateEvaluateCovariance(Covaₗₘ)
    EvaluateFisherMatrix!(VariedParameters, Fisher, Path∂Cℓ, CovCℓ)
    SelectMatrixAndMarginalize!(VariedParameters, Fisher)
    return Fisher
end