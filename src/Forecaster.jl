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
                w0waCDMCosmology = CosmoCentral.w0waCDMCosmology(
                w0 = CopyDictCosmo["w0"][1], wa = CopyDictCosmo["wa"][1],
                Mν = CopyDictCosmo["Mν"][1], H0 = CopyDictCosmo["H0"][1],
                ΩM = CopyDictCosmo["ΩM"][1], ΩB = CopyDictCosmo["ΩB"][1],
                ΩDE = CopyDictCosmo["ΩDE"][1], Ωk = CopyDictCosmo["Ωk"][1],
                Ωr = CopyDictCosmo["Ωr"][1], ns = CopyDictCosmo["ns"][1],
                σ8 = CopyDictCosmo["σ8"][1])
                MyDict["dvar_"*key*"_step_p_"*string(index)] = [w0waCDMCosmology]
                myvalue = IncrementedValue(value[1], -mystep)
                CopyDictCosmo[key] = [myvalue]
                w0waCDMCosmology = CosmoCentral.w0waCDMCosmology(w0 =
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
    w0waCDMCosmology = CosmoCentral.w0waCDMCosmology(
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
    LensingFunction::WLWeightFunction)
    ComputeLensingEfficiencyGrid!(LensingFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology, CustomLensingEfficiency())
    ComputeWeightFunctionGrid!(LensingFunction, ConvolvedDensity,
    CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
    return LensingFunction
end

function Forecast∂Cℓ!(DictCosmo::Dict,
    PathInput::String, PathConfig::String, Steps::Array)
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
    end
end

function InstantiateWL(DictInput::Dict)
    LensingFunction = WLWeightFunction()
    return LensingFunction
end

function InstantiateGC(DictInput::Dict)
    gcWeightFunction = GCWeightFunction()
    InstantiateBias(DictInput, gcWeightFunction)
    return gcWeightFunction
end

function InstantiateBias(DictInput::Dict,
    GCWeightFunction::GCWeightFunction)
    if DictInput["PhotometricGalaxy"]["Bias"] == "PiecewiseBias"
        GCWeightFunction.BiasKind = PiecewiseBias()
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
        LensingFunction)
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
            #TODO is this if-else necessary?
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

function ForecastCℓ!(Cosmologies::Dict, PathInput::String,
    PathOutput::String, CosmologicalGrid::CosmologicalGrid, PathConfig::String)
    ProbesDict = JSON.parsefile(PathConfig)
    analyticaldensity = AnalitycalDensity()
    NormalizeAnalitycalDensity!(analyticaldensity)
    instrumentresponse = InstrumentResponse()
    convolveddensity = ConvolvedDensity(DensityGridArray =
    ones(10, length(CosmologicalGrid.ZArray)))
    NormalizeConvolvedDensity!(convolveddensity, analyticaldensity,
    instrumentresponse, CosmologicalGrid)
    ComputeConvolvedDensityGrid!(CosmologicalGrid, convolveddensity,
    analyticaldensity, instrumentresponse)
    for (key, value) in Cosmologies
        w0waCDMCosmology = value[1]
        PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
        ReadPowerSpectrumBackground(PathInput*key*"/p_mm",
        CosmologicalGrid.MultipolesArray)
        CosmologicalGrid.MultipolesArray[1,1]
        DictProbes = InitializeProbes(ProbesDict, convolveddensity,
        w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities)
        ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
        InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
        BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
        InitializeForecastCℓ(DictProbes, BackgroundQuantities,
        w0waCDMCosmology, CosmologicalGrid, PowerSpectrum, PathOutput, key)
    end
end
