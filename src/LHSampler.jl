function CosmologicalParameterSampler(DictCosmo::Dict, n_points::Int)
    parameter_list = ["w0", "wa", "Mν", "H0", "ΩM", "ΩB", "ΩDE", "Ωk", "Ωr", "ns", "σ8"]
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
        w0waCDMCosmology = CosmoCentral.w0waCDMCosmology(
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
        BackgroundQuantities = BackgroundQuantities(
        HZArray = zeros(length(CosmologicalGrid.ZArray)),
        rZArray=zeros(length(CosmologicalGrid.ZArray)))
        ComputeBackgroundQuantitiesGrid!(CosmologicalGrid,
        BackgroundQuantities, value[1])
        ClassyParams = Initializeclassy(value[1])
        PowerSpectrum = PowerSpectrum(PowerSpectrumLinArray =
        zeros(length(CosmologicalGrid.KArray), length(CosmologicalGrid.ZArray)),
        PowerSpectrumNonlinArray = zeros(length(CosmologicalGrid.KArray),
        length(CosmologicalGrid.ZArray)),
        InterpolatedPowerSpectrum = zeros(length(
        CosmologicalGrid.ℓBinCenters), length(CosmologicalGrid.ZArray)))
        EvaluatePowerSpectrum!(ClassyParams, CosmologicalGrid, PowerSpectrum)
        WritePowerSpectrumBackground(PowerSpectrum, BackgroundQuantities,
        CosmologicalGrid, Path*key*"/p_mm")
        WriteCosmology!(value[1], Path*key)
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

function EvaluateCℓLHS!(Cosmologies::Dict, PathInput::String,
    PathOutput::String, CosmologicalGrid::CosmologicalGrid, PathConfig::String)
    ProbesDict = JSON.parsefile(PathConfig)
    AnalitycalDensity = AnalitycalDensity()
    NormalizeAnalitycalDensity!(AnalitycalDensity)
    InstrumentResponse = InstrumentResponse()
    ConvolvedDensity = ConvolvedDensity(DensityGridArray =
    ones(10, length(CosmologicalGrid.ZArray)))
    NormalizeConvolvedDensity!(ConvolvedDensity, AnalitycalDensity,
    InstrumentResponse, CosmologicalGrid)
    ComputeConvolvedDensityGrid!(CosmologicalGrid, ConvolvedDensity,
    AnalitycalDensity, InstrumentResponse)
    for (key, value) in Cosmologies
        w0waCDMCosmology = value[1]
        PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
        ReadPowerSpectrumBackground(PathInput*key*"/p_mm",
        CosmologicalGrid.ℓBinCenters)
        CosmologicalGrid.ℓBinCenters[1,1]
        DictProbes = InitializeProbes(ProbesDict, ConvolvedDensity,
        w0waCDMCosmology, CosmologicalGrid, BackgroundQuantities)
        ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
        InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
        BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
        InitializeForecastCℓ(DictProbes, BackgroundQuantities,
        w0waCDMCosmology, CosmologicalGrid, PowerSpectrum, PathOutput, key)
    end
end

function EvaluateCℓGeneral!(PmmDirectory::String,
    PathOutput::String, CosmologicalGrid::CosmologicalGrid, PathConfig::String,
    PathInput::String)
    ProbesDict = JSON.parsefile(PathConfig)
    analitycaldensity = AnalitycalDensity()
    NormalizeAnalitycalDensity!(analitycaldensity)
    instrumentresponse = InstrumentResponse()
    convolveddensity = ConvolvedDensity(DensityGridArray =
    ones(10, length(CosmologicalGrid.ZArray)))
    NormalizeConvolvedDensity!(convolveddensity, analitycaldensity,
    instrumentresponse, CosmologicalGrid)
    ComputeConvolvedDensityGrid!(CosmologicalGrid, convolveddensity,
    analitycaldensity, instrumentresponse)
    for (root, dirs, files) in walkdir(PmmDirectory)
        for file in files
            file_extension = file[findlast(isequal('.'),file):end]
            if file_extension == ".json"
                CosmoDict = JSON.parsefile(joinpath(root, file))
                w0waCDMCosmology = CosmoCentral.ReadCosmology(Dict(CosmoDict))
                PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
                ReadPowerSpectrumBackground(joinpath(root, "p_mm"),
                CosmologicalGrid.ℓBinCenters)
                ExtractGrowthFactor!(BackgroundQuantities, PowerSpectrum)
                #TODO probably we can obtain flexibility with a dictionary
                intrinsicalignment = ReadIntrinsicAlignment(Dict(CosmoDict))
                bias = ReadBias(Dict(CosmoDict))
                CopyConvolvedDensity = deepcopy(convolveddensity)
                #CopyConvolvedDensity.ShiftArray = ones(10).*CosmoDict["ShiftParameter"]
                #ShiftConvolvedDensityGrid!(CosmologicalGrid, CopyConvolvedDensity)
                DictProbes = InitializeProbes(ProbesDict, CopyConvolvedDensity,
                w0waCDMCosmology, intrinsicalignment, bias, CosmologicalGrid,
                BackgroundQuantities, PathInput)
                ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
                InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
                BackgroundQuantities, PowerSpectrum, BSplineCubic())
                RandomString = Random.randstring(12)
                mkdir(joinpath(PathOutput,RandomString))
                InitializeForecastCℓ(DictProbes,
                BackgroundQuantities, w0waCDMCosmology, CosmologicalGrid,
                PowerSpectrum, PathOutput, CosmoDict, RandomString)
            end
        end
    end
end

function EvaluateCℓDoubleShift(PmmDirectory::String,
    PathOutput::String, CosmologicalGrid::CosmologicalGrid, PathConfig::String)
    ProbesDict = JSON.parsefile(PathConfig)
    AnalitycalDensity = AnalitycalDensity()
    NormalizeAnalitycalDensity!(AnalitycalDensity)
    InstrumentResponse = InstrumentResponse()
    ConvolvedDensity = ConvolvedDensity(DensityGridArray =
    ones(10, length(CosmologicalGrid.ZArray)))
    NormalizeConvolvedDensity!(ConvolvedDensity, AnalitycalDensity,
    InstrumentResponse, CosmologicalGrid)
    ComputeConvolvedDensityGrid!(CosmologicalGrid, ConvolvedDensity,
    AnalitycalDensity, InstrumentResponse)
    #TODO Probably the following line is useless
    ℓBinCenters = Array(LinRange(10,3000,100))
    for (root, dirs, files) in walkdir(PmmDirectory)
        for file in files
            file_extension = file[findlast(isequal('.'),file):end]
            if file_extension == ".json"
                CosmoDict = JSON.parsefile(joinpath(root, file))
                w0waCDMCosmology = CosmoCentral.ReadCosmology(Dict(CosmoDict))
                PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
                ReadPowerSpectrumBackground(joinpath(root, "p_mm"),
                CosmologicalGrid.ℓBinCenters)
                ShiftParameter_i = CosmoDict["ShiftParameter_i"]
                ShiftParameter_j = CosmoDict["ShiftParameter_j"]
                ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
                InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
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
                            cℓ = Cℓ(
                            CℓArray
                            = zeros(length(CosmologicalGrid.ℓBinCenters),
                            10, 10))
                            push!(CℓArray, cℓ)
                        end
                    end
                end
                for i in 1:10
                    for j in i+1:10
                        CopyConvolvedDensity = deepcopy(ConvolvedDensity)
                        CopyConvolvedDensity.ShiftArray[i] = ShiftParameter_i
                        CopyConvolvedDensity.ShiftArray[j] = ShiftParameter_j
                        ShiftConvolvedDensityGrid!(CosmologicalGrid,
                        CopyConvolvedDensity)
                        DictProbes = InitializeProbes(ProbesDict,
                        CopyConvolvedDensity, w0waCDMCosmology,
                        CosmologicalGrid, BackgroundQuantities)
                        CℓKeyArray, TempCℓArray =
                        InitializeComputeCℓDoubleNuisance(
                        DictProbes, BackgroundQuantities, w0waCDMCosmology,
                        CosmologicalGrid, PowerSpectrum)
                        for (index, myCℓ ) in enumerate(TempCℓArray)
                            (CℓArray[index]).CℓArray[:,i,j] .=
                            myCℓ.CℓArray[:,i,j]
                            (CℓArray[index]).CℓArray[:,j,i] .=
                            myCℓ.CℓArray[:,j,i]
                        end
                    end
                end
                RandomString = Random.randstring(12)
                mkdir(joinpath(PathOutput,RandomString))
                WriteCℓ!("Lensing_Lensing",
                CℓArray[1], PathOutput*RandomString*"/cl")
                WriteCℓ!("Lensing_PhotometricGalaxy",
                CℓArray[2], PathOutput*RandomString*"/cl")
                WriteCℓ!("PhotometricGalaxy_PhotometricGalaxy",
                CℓArray[3], PathOutput*RandomString*"/cl")
                #TODO this function is not beautiful. It is quite long and it
                # needs a serious refactoring
                WriteParameters!(CosmoDict, PathOutput*RandomString)
            end
        end
    end

end

function InitializeComputeCℓDoubleNuisance(ProbesDict::Dict,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology,
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
                cℓ = Cℓ(
                CℓArray
                = zeros(length(CosmologicalGrid.ℓBinCenters),
                length(ProbesDict[key_A].WeightFunctionArray[:, 1]),
                length(ProbesDict[key_B].WeightFunctionArray[:, 1])))
                ComputeCℓ!(cℓ,
                ProbesDict[key_A], ProbesDict[key_B], BackgroundQuantities,
                w0waCDMCosmology, CosmologicalGrid, PowerSpectrum,
                CosmoCentral.CustomSimpson())
                push!(CℓArray, cℓ)
            end
        end
    end
    return CℓKeyArray, CℓArray
end
