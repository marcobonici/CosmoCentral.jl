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

function EvaluateAngularCoefficientsGeneral(PmmDirectory::String,
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
    #TODO Probably the following line is useless
    MultipolesArray = Array(LinRange(10,3000,100))
    for (root, dirs, files) in walkdir(PmmDirectory)
        for file in files
            file_extension = file[findlast(isequal('.'),file):end]
            if file_extension == ".json"
                CosmoDict = JSON.parsefile(joinpath(root, file))
                w0waCDMCosmology = CosmoCentral.ReadCosmology(Dict(CosmoDict))
                PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
                ReadPowerSpectrumBackground(joinpath(root, "p_mm"),
                CosmologicalGrid.MultipolesArray)
                CopyConvolvedDensity = deepcopy(ConvolvedDensity)
                CopyConvolvedDensity.ShiftArray =
                ones(10).*CosmoDict["ShiftParameter"]
                ShiftConvolvedDensityFunctionGrid(CosmologicalGrid,
                CopyConvolvedDensity)
                DictProbes = InitializeProbes(ProbesDict, CopyConvolvedDensity,
                w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities)
                ComputeLimberArray(CosmologicalGrid, BackgroundQuantities)
                InterpolateAndEvaluatePowerSpectrum(CosmologicalGrid,
                BackgroundQuantities, PowerSpectrum, BSplineCubic())
                RandomString = Random.randstring(12)
                mkdir(joinpath(PathOutput,RandomString))
                InitializeComputeAngularCoefficients(DictProbes,
                BackgroundQuantities, w0waCDMCosmology, CosmologicalGrid,
                PowerSpectrum, PathOutput, CosmoDict, RandomString)
            end
        end
    end
end

function EvaluateAngularCoefficientsDoubleShift(PmmDirectory::String,
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
    #TODO Probably the following line is useless
    MultipolesArray = Array(LinRange(10,3000,100))
    for (root, dirs, files) in walkdir(PmmDirectory)
        for file in files
            file_extension = file[findlast(isequal('.'),file):end]
            if file_extension == ".json"
                CosmoDict = JSON.parsefile(joinpath(root, file))
                w0waCDMCosmology = CosmoCentral.ReadCosmology(Dict(CosmoDict))
                PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
                ReadPowerSpectrumBackground(joinpath(root, "p_mm"),
                CosmologicalGrid.MultipolesArray)
                ShiftParameter_i = CosmoDict["ShiftParameter_i"]
                ShiftParameter_j = CosmoDict["ShiftParameter_j"]
                ComputeLimberArray(CosmologicalGrid, BackgroundQuantities)
                InterpolateAndEvaluatePowerSpectrum(CosmologicalGrid,
                BackgroundQuantities, PowerSpectrum, BSplineCubic())
                ProbesArray = []
                CℓKeyArray = []
                CℓArray = []
                for (key, value) in ProbesDict
                    push!(ProbesArray, key)
                end
                sort!(ProbesArray)
                for key_A in ProbesArray
                    for key_B in ProbesArray
                        if key_B*"_"*key_A in CℓKeyArray
                        else
                            push!(CℓKeyArray, key_A*"_"*key_B)
                            AngularCoefficients = AngularCoefficientsStruct(
                            AngularCoefficientsArray
                            = zeros(length(CosmologicalGrid.MultipolesArray),
                            10, 10))
                            push!(CℓArray, AngularCoefficients)
                        end
                    end
                end
                for i in 1:10
                    for j in i+1:10
                        CopyConvolvedDensity = deepcopy(ConvolvedDensity)
                        CopyConvolvedDensity.ShiftArray[i] = ShiftParameter_i
                        CopyConvolvedDensity.ShiftArray[j] = ShiftParameter_j
                        ShiftConvolvedDensityFunctionGrid(CosmologicalGrid,
                        CopyConvolvedDensity)
                        DictProbes = InitializeProbes(ProbesDict,
                        CopyConvolvedDensity, w0waCDMCosmology,
                        CosmologicalGrid, BackgroundQuantities)
                        CℓKeyArray, TempCℓArray =
                        InitializeComputeAngularCoefficientsDoubleNuisance(
                        DictProbes, BackgroundQuantities, w0waCDMCosmology,
                        CosmologicalGrid, PowerSpectrum)
                        for (index, myCℓ ) in enumerate(TempCℓArray)
                            (CℓArray[index]).AngularCoefficientsArray[:,i,j] .=
                            myCℓ.AngularCoefficientsArray[:,i,j]
                            (CℓArray[index]).AngularCoefficientsArray[:,j,i] .=
                            myCℓ.AngularCoefficientsArray[:,j,i]
                        end
                    end
                end
                RandomString = Random.randstring(12)
                mkdir(joinpath(PathOutput,RandomString))
                WriteAngularCoefficients(key_A*"_"*key_B,
                AngularCoefficients, PathOutput*RandomString*"/cl")
                WriteParameters(CosmoDict, PathOutput*RandomString)
            end
        end
    end

end

function InitializeComputeAngularCoefficientsDoubleNuisance(ProbesDict::Dict,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    CosmologicalGrid::CosmologicalGrid, PowerSpectrum::PowerSpectrum)
    ProbesArray = []
    CℓKeyArray = []
    CℓArray = []
    for (key, value) in ProbesDict
        push!(ProbesArray, key)
    end
    sort!(ProbesArray)
    for key_A in ProbesArray
        for key_B in ProbesArray
            if key_B*"_"*key_A in CℓKeyArray
            else
                push!(CℓKeyArray, key_A*"_"*key_B)
                AngularCoefficients = AngularCoefficientsStruct(
                AngularCoefficientsArray
                = zeros(length(CosmologicalGrid.MultipolesArray),
                length(ProbesDict[key_A].WeightFunctionArray[:, 1]),
                length(ProbesDict[key_B].WeightFunctionArray[:, 1])))
                ComputeAngularCoefficients(AngularCoefficients,
                ProbesDict[key_A], ProbesDict[key_B], BackgroundQuantities,
                w0waCDMCosmology, CosmologicalGrid, PowerSpectrum,
                CosmoCentral.CustomTrapz())
                push!(CℓArray, AngularCoefficients)
            end
        end
    end
    return CℓKeyArray, CℓArray
end
