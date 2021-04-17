"""
    ComputeWeightFunction(z::Float64, i::Int64,
        ConvolvedDensity::ConvolvedDensity,
        w0waCDMCosmology::w0waCDMCosmologyStruct)

This function returns the Galaxy Clustering Weight function for a given redshift
``z`` and tomographic bin ``i``.
"""
function ComputeWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::AsbtractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    GCWeightFunction::GCWeightFunctionStruct)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    return ComputeConvolvedDensityFunction(z, i, ConvolvedDensity,
    AnalitycalDensity, InstrumentResponse) *
    ComputeBias(z, GCWeightFunction.BiasKind, ConvolvedDensity) *
    ComputeHubbleFactor(z, w0waCDMCosmology) / c_0
end

"""
    ComputeWeightFunctionOverGrid(GCWeightFunction::GCWeightFunctionStruct,
    ConvolvedDensity::AsbtractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct)

This function evaluates the Galaxy Clustering Weight Function over the ``z``
grid and for all tomographic bins ``i``.
"""
function ComputeWeightFunctionOverGrid(GCWeightFunction::GCWeightFunctionStruct,
    ConvolvedDensity::AsbtractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    #TODO use avx here to accelerate
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            GCWeightFunction.WeightFunctionArray[idx_ZBinArray, idx_ZArray] =
            GCWeightFunction.BiasArray[idx_ZBinArray, idx_ZArray]*
            ConvolvedDensity.DensityGridArray[idx_ZBinArray,
            idx_ZArray] *
            BackgroundQuantities.HZArray[idx_ZArray] / c_0
        end
    end
end

"""
    ComputeLensingEfficiency(z::Float64, i::Int64,
    ConvolvedDensity::AsbtractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    CosmologicalGrid::CosmologicalGrid,
    WLWeightFunction::WLWeightFunctionStruct)

This function returns the Lensing efficiency, for a given redshift
``z`` and tomographic bin ``i``.
"""
function ComputeLensingEfficiency(z::Float64, i::Int64,
    ConvolvedDensity::AsbtractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    CosmologicalGrid::CosmologicalGrid,
    WLWeightFunction::WLWeightFunctionStruct)
    int, err = QuadGK.quadgk(x -> CosmoCentral.ComputeConvolvedDensityFunction(
    x, i, ConvolvedDensity, AnalitycalDensity, InstrumentResponse)*
    (1. - CosmoCentral.ComputeComovingDistance(z, w0waCDMCosmology)/
    CosmoCentral.ComputeComovingDistance(x, w0waCDMCosmology)), z,
    last(CosmologicalGrid.ZArray) , rtol=1e-12)
    return int
end

function ComputeLensingEfficiency(z::Float64, i::Int64,
    ConvolvedDensity::AsbtractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    CosmologicalGrid::CosmologicalGrid,
    LensingSourceFunction::LensingSourceFunctionStruct)
    int, err = QuadGK.quadgk(x -> CosmoCentral.ComputeConvolvedDensityFunction(
    x, i, ConvolvedDensity, AnalitycalDensity, InstrumentResponse)*
    ((CosmoCentral.ComputeComovingDistance(x, w0waCDMCosmology) - CosmoCentral.ComputeComovingDistance(z, w0waCDMCosmology))/
    (CosmoCentral.ComputeComovingDistance(x, w0waCDMCosmology) * CosmoCentral.ComputeComovingDistance(z, w0waCDMCosmology)  ) ), z,
    last(CosmologicalGrid.ZArray) , rtol=1e-12)
    return int
end


"""
    ComputeLensingEfficiencyOverGrid(
    WLWeightFunction::WLWeightFunctionStruct,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse,
    ConvolvedDensity::AsbtractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct)

This function evaluates the Lensing Efficiency over the ``z``
grid and for all tomographic bins ``i``.
"""
function ComputeLensingEfficiencyOverGrid(
    WLWeightFunction::WLWeightFunctionStruct,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse,
    ConvolvedDensity::AsbtractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct)
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            WLWeightFunction.LensingEfficiencyArray[idx_ZBinArray, idx_ZArray] =
            ComputeLensingEfficiency(CosmologicalGrid.ZArray[idx_ZArray],
            idx_ZBinArray, ConvolvedDensity, AnalitycalDensity,
            InstrumentResponse, w0waCDMCosmology, CosmologicalGrid,
            WLWeightFunction)
        end
    end
end

function ComputeLensingEfficiencyOverGrid(
    LensingSourceFunction::LensingSourceFunctionStruct,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse,
    ConvolvedDensity::AsbtractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct)
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            LensingSourceFunction.LensingEfficiencyArray[idx_ZBinArray, idx_ZArray] =
            ComputeLensingEfficiency(CosmologicalGrid.ZArray[idx_ZArray],
            idx_ZBinArray, ConvolvedDensity, AnalitycalDensity,
            InstrumentResponse, w0waCDMCosmology, CosmologicalGrid,
            LensingSourceFunction)
        end
    end
end

"""
    ComputeLensingEfficiencyOverGridCustom(
    WLWeightFunction::WLWeightFunctionStruct,
    ConvolvedDensity::AsbtractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct)

This function evaluates the Lensing Efficiency over the ``z``
grid and for all tomographic bins ``i``.
"""
function ComputeLensingEfficiencyOverGridCustom(
    WLWeightFunction::WLWeightFunctionStruct,
    ConvolvedDensity::AsbtractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct)
    WLWeightFunction.LensingEfficiencyArray .*= 0
    Weight_Matrix = SimpsonWeightMatrix(length(CosmologicalGrid.ZArray))
    @avx for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            for idx_ZArrayInt in 1:length(CosmologicalGrid.ZArray)
                WLWeightFunction.LensingEfficiencyArray[idx_ZBinArray,
                idx_ZArray] += ConvolvedDensity.DensityGridArray[idx_ZBinArray,
                idx_ZArrayInt] * (1 - BackgroundQuantities.rZArray[idx_ZArray] /
                BackgroundQuantities.rZArray[idx_ZArrayInt]) *
                Weight_Matrix[idx_ZArray, idx_ZArrayInt]
            end
        end
    end
    WLWeightFunction.LensingEfficiencyArray .*=
    (last(CosmologicalGrid.ZArray)-first(CosmologicalGrid.ZArray))/
    (length(CosmologicalGrid.ZArray)-1)
end

"""
    ComputeWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::AsbtractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    CosmologicalGrid::CosmologicalGrid,
    WLWeightFunction::WLWeightFunctionStruct)

This function returns the Weak Lensing Weight Function, for a given redshift
``z`` and tomographic bin ``i``.
"""
function ComputeWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::AsbtractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse,
    w0waCDMCosmology::w0waCDMCosmologyStruct,
    CosmologicalGrid::CosmologicalGrid,
    WLWeightFunction::WLWeightFunctionStruct)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    return 1.5 * ComputeLensingEfficiency(z, i, ConvolvedDensity,
    AnalitycalDensity,InstrumentResponse, w0waCDMCosmology, CosmologicalGrid,
    WLWeightFunction) * (w0waCDMCosmology.H0/c_0)^2 * w0waCDMCosmology.ΩM *
    (1. + z) * ComputeComovingDistance(z, w0waCDMCosmology)
end


"""
    ComputeWeightFunctionOverGrid(
    WLWeightFunction::WLWeightFunctionStruct,
    ConvolvedDensity::AsbtractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct)

This function evaluates the Weak Lensing Weight Function over the ``z``
grid and for all tomographic bins ``i``.
"""
function ComputeWeightFunctionOverGrid(
    WLWeightFunction::WLWeightFunctionStruct,
    ConvolvedDensity::AsbtractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmologyStruct)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            WLWeightFunction.WeightFunctionArray[idx_ZBinArray, idx_ZArray] =
            1.5 * (w0waCDMCosmology.H0/c_0)^2 * w0waCDMCosmology.ΩM *
            (1. + CosmologicalGrid.ZArray[idx_ZArray]) *
            BackgroundQuantities.rZArray[idx_ZArray] *
            WLWeightFunction.LensingEfficiencyArray[idx_ZBinArray, idx_ZArray]
        end
    end
end
