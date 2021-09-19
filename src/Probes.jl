abstract type AbstractProbe end

@kwdef mutable struct GCProbe <: AbstractProbe
    Density::AbstractConvolvedDensity
    BiasModel::AbstractBias
    NuisanceDict::Dict = Dict()
end

@kwdef mutable struct WLProbe <:AbstractProbe
    Density::AbstractConvolvedDensity
    IAModel::AbstractIntrinsicAlignment
    NuisanceDict::Dict = Dict()
end

function CreateWLProbe(DensityDict::Dict, IADict::Dict, cosmogrid::CosmologicalGrid)
    wlprobe = WLProbe(Density = CreateDensity(DensityDict, cosmogrid),
    IAModel = CreateIA(IADict))
    return wlprobe
end

function CreateGCProbe(DensityDict::Dict, BiasDict::Dict, cosmogrid::CosmologicalGrid)
    gcprobe = GCProbe(Density = CreateDensity(DensityDict, cosmogrid),
    BiasModel = CreateBias(BiasDict))
    return gcprobe
end

function IterateWLProbe(ProbesDict::Dict, cosmogrid::CosmoCentral.CosmologicalGrid,
    SteMSteps::Array)
    LensingProbeDict = Dict()
    IADict = ReadPresentDict(ProbesDict["IntrinsicAlignment"])
    LensingProbeDict["dvar_central_step_0"] = CreateWLProbe(ProbesDict["density"], IADict, cosmogrid)
    TempIADict = Dict()
    Iterator!(TempIADict, ProbesDict["IntrinsicAlignment"], SteMSteps)
    for (key,value) in TempIADict
        LensingProbeDict[key] = CreateWLProbe(ProbesDict["density"], value, cosmogrid)
    end
    return LensingProbeDict
end

function IterateGCProbe(ProbesDict::Dict, cosmogrid::CosmoCentral.CosmologicalGrid,
    SteMSteps::Array)
    GCProbeDict = Dict()
    BiasDict = ReadPresentDict(ProbesDict["Bias"])
    GCProbeDict["dvar_central_step_0"] = CreateGCProbe(ProbesDict["density"], BiasDict, cosmogrid)
    TempBiasDict = Dict()
    Iterator!(TempBiasDict, ProbesDict["Bias"], SteMSteps)
    for (key,value) in TempBiasDict
        GCProbeDict[key] = CreateGCProbe(ProbesDict["density"], value, cosmogrid)
    end
    return GCProbeDict
end

function IterateProbe(ProbesDict::Dict, cosmogrid::CosmoCentral.CosmologicalGrid, SteMSteps::Array)
    if ProbesDict["probe"] == "WLProbe"
        IteratedDict = IterateWLProbe(ProbesDict, cosmogrid, SteMSteps)
    elseif ProbesDict["probe"] == "GCProbe"
        IteratedDict = IterateGCProbe(ProbesDict, cosmogrid, SteMSteps)
    else
        error("No probes!")
    end
    return IteratedDict
end
