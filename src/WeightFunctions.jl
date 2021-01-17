"""
    ComputeWeightFunction(z::Float64, i::Int64,
        ConvolvedDensity::ConvolvedDensity,
        w0waCDMCosmology::w0waCDMCosmology)

This function returns the Galaxy Clustering Weight function for a given redshift
``z`` and tomographic bin ``i``.
"""
function ComputeWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::ConvolvedDensity, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, w0waCDMCosmology::w0waCDMCosmology,
    Bias::Bias, GCWeightFunction::GCWeightFunction)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    return ComputeConvolvedDensityFunction(z, i, ConvolvedDensity,
    AnalitycalDensity, InstrumentResponse) *
    ComputeBias(z, Bias, ConvolvedDensity) *
    ComputeHubbleFactor(z, w0waCDMCosmology) / c_0
end

"""
    ComputeWeightFunctionOverGrid(
    GCWeightFunction::GCWeightFunction, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, ConvolvedDensity::ConvolvedDensity,
    Bias::Bias, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology)

This function evaluates the Galaxy Clustering Weight Function over the ``z``
grid and for all tomographic bins ``i``.
"""
function ComputeWeightFunctionOverGrid(
    GCWeightFunction::GCWeightFunction, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, ConvolvedDensity::ConvolvedDensity,
    Bias::Bias, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            GCWeightFunction.WeightFunctionArray[idx_ZBinArray, idx_ZArray] =
            Bias.BiasArray[idx_ZBinArray, idx_ZArray]*
            ConvolvedDensity.DensityGridArray[idx_ZBinArray,
            idx_ZArray] *
            BackgroundQuantities.HZArray[idx_ZArray] / c_0
        end
    end
end

"""
    ComputeWeightFunction(z::Float64, i::Int64,
        ConvolvedDensity::ConvolvedDensity,
        w0waCDMCosmology::w0waCDMCosmology)

This function returns the Galaxy Clustering Weight function for a given redshift
``z`` and tomographic bin ``i``.
"""
function ComputeLensingEfficiency(z::Float64, i::Int64,
    ConvolvedDensity::ConvolvedDensity, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, w0waCDMCosmology::w0waCDMCosmology,
    CosmologicalGrid::CosmologicalGrid, WLWeightFunction::WLWeightFunction)
    int, err = QuadGK.quadgk(x -> CosmoCentral.ComputeConvolvedDensityFunction(x, i,
    ConvolvedDensity, AnalitycalDensity, InstrumentResponse)*
    (1. - CosmoCentral.ComputeComovingDistance(z, w0waCDMCosmology)/
    CosmoCentral.ComputeComovingDistance(x, w0waCDMCosmology)), z,
    last(CosmologicalGrid.ZArray) , rtol=1e-12)
    return int
end

"""
    ComputeWeightFunctionOverGrid(
    GCWeightFunction::GCWeightFunction, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, ConvolvedDensity::ConvolvedDensity,
    Bias::Bias, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology)

This function evaluates the Galaxy Clustering Weight Function over the ``z``
grid and for all tomographic bins ``i``.
"""
function ComputeLensingEfficiencyOverGrid(
    WLWeightFunction::WLWeightFunction, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, ConvolvedDensity::ConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology)
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

function ComputeWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::ConvolvedDensity, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, w0waCDMCosmology::w0waCDMCosmology,
    CosmologicalGrid::CosmologicalGrid, WLWeightFunction::WLWeightFunction)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    return 1.5 * ComputeLensingEfficiency(z, i, ConvolvedDensity,
    AnalitycalDensity,InstrumentResponse, w0waCDMCosmology, CosmologicalGrid,
    WLWeightFunction) * (w0waCDMCosmology.H0/c_0)^2 * w0waCDMCosmology.ΩM *
    (1. + z) * ComputeComovingDistance(z, w0waCDMCosmology)
end


function ComputeWeightFunctionOverGrid(
    WLWeightFunction::WLWeightFunction, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, ConvolvedDensity::ConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology)
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
