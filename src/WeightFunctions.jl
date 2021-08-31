"""
    ComputeWeightFunction(z::Float64, i::Int64, ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity, InstrumentResponse::InstrumentResponse,
    w0waCDMCosmology::w0waCDMCosmology, GCWeightFunction::GCWeightFunction)

This function returns the Galaxy Clustering Weight function for a given redshift
``z`` and tomographic bin ``i``.
"""
function ComputeWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::AbstractConvolvedDensity, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, Cosmology::AbstractCosmology,
    GCWeightFunction::GCWeightFunction)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    return ComputeConvolvedDensity(z, i, ConvolvedDensity,
    AnalitycalDensity, InstrumentResponse) *
    ComputeBias(z, GCWeightFunction.BiasKind, ConvolvedDensity) *
    ComputeHubbleFactor(z, Cosmology) / c_0
end

"""
    ComputeWeightFunctionGrid!(GCWeightFunction::GCWeightFunction,
    ConvolvedDensity::AbstractConvolvedDensity, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, Cosmology::AbstractCosmology)

This function evaluates the Galaxy Clustering Weight Function over the ``z``
grid and for all tomographic bins ``i``.
"""
function ComputeWeightFunctionGrid!(GCWeightFunction::GCWeightFunction,
    ConvolvedDensity::AbstractConvolvedDensity, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, Cosmology::AbstractCosmology)
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
    ConvolvedDensity::AbstractConvolvedDensity, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, w0waCDMCosmology::w0waCDMCosmology,
    CosmologicalGrid::CosmologicalGrid, ::WLWeightFunction)

This function returns the Lensing efficiency, for a given redshift
``z`` and tomographic bin ``i``.
"""
function ComputeLensingEfficiency(z::Float64, i::Int64,
    ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse,
    Cosmology::AbstractCosmology,
    CosmologicalGrid::CosmologicalGrid,
    ::WLWeightFunction)
    int, err = QuadGK.quadgk(x -> CosmoCentral.ComputeConvolvedDensity(
    x, i, ConvolvedDensity, AnalitycalDensity, InstrumentResponse)*
    ((Computeχ(x, Cosmology) -
    Computeχ(z, Cosmology))/
    (Computeχ(x, Cosmology) *
    Computeχ(z, Cosmology)  ) ), z,
    last(CosmologicalGrid.ZArray) , rtol=1e-12)
    return int
end

function ComputeLensingEfficiency(z::Float64, i::Int64,
    ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse,
    Cosmology::AbstractCosmology,
    CosmologicalGrid::CosmologicalGrid,
    LensingSourceFunction::LensingSourceFunction)
    int, err = QuadGK.quadgk(x -> CosmoCentral.ComputeConvolvedDensity(
    x, i, ConvolvedDensity, AnalitycalDensity, InstrumentResponse)*
    ((Computeχ(x, Cosmology) -
    Computeχ(z, Cosmology))/
    (Computeχ(x, Cosmology) *
    Computeχ(z, Cosmology)  ) ), z,
    last(CosmologicalGrid.ZArray) , rtol=1e-12)
    return int
end


"""
    ComputeLensingEfficiencyGrid!(LensingFunction::WLWeightFunction,
    AnalitycalDensity::AnalitycalDensity, InstrumentResponse::InstrumentResponse,
    ConvolvedDensity::AbstractConvolvedDensity, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, Cosmology::AbstractCosmology,
    ::StandardLensingEfficiency)

This function evaluates the Lensing Efficiency over the ``z``
grid and for all tomographic bins ``i``.
"""
function ComputeLensingEfficiencyGrid!(LensingFunction::WLWeightFunction,
    AnalitycalDensity::AnalitycalDensity, InstrumentResponse::InstrumentResponse,
    ConvolvedDensity::AbstractConvolvedDensity, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, Cosmology::AbstractCosmology,
    ::StandardLensingEfficiency)
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            LensingFunction.LensingEfficiencyArray[idx_ZBinArray, idx_ZArray] =
            ComputeLensingEfficiency(CosmologicalGrid.ZArray[idx_ZArray],
            idx_ZBinArray, ConvolvedDensity, AnalitycalDensity,
            InstrumentResponse, Cosmology, CosmologicalGrid,
            LensingFunction)
        end
    end
end

function ComputeLensingEfficiencyGrid!(
    LensingSourceFunction::LensingSourceFunction,
    AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse,
    ConvolvedDensity::AbstractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    Cosmology::AbstractCosmology,
    ::StandardLensingEfficiency)
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            LensingSourceFunction.LensingEfficiencyArray[idx_ZBinArray, idx_ZArray] =
            ComputeLensingEfficiency(CosmologicalGrid.ZArray[idx_ZArray],
            idx_ZBinArray, ConvolvedDensity, AnalitycalDensity,
            InstrumentResponse, Cosmology, CosmologicalGrid,
            LensingSourceFunction)
        end
    end
end

"""
    ComputeLensingEfficiencyGrid!(
    LensingFunction::WLWeightFunction, ConvolvedDensity::AbstractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid, BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology, ::CustomLensingEfficiency)

This function evaluates the Lensing Efficiency over the ``z`` grid and for all tomographic
bins ``i``.
"""
function ComputeLensingEfficiencyGrid!(
    LensingFunction::WLWeightFunction,
    ConvolvedDensity::AbstractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    Cosmology::AbstractCosmology,
    ::CustomLensingEfficiency)
    LensingFunction.LensingEfficiencyArray .*= 0
    Weight_Matrix = SimpsonWeightMatrix(length(CosmologicalGrid.ZArray))
    @avx for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            for idx_ZArrayInt in 1:length(CosmologicalGrid.ZArray)
                LensingFunction.LensingEfficiencyArray[idx_ZBinArray,
                idx_ZArray] += ConvolvedDensity.DensityGridArray[idx_ZBinArray,
                idx_ZArrayInt] * (BackgroundQuantities.χZArray[idx_ZArrayInt] -
                BackgroundQuantities.χZArray[idx_ZArray]) /
                BackgroundQuantities.χZArray[idx_ZArrayInt] /
                BackgroundQuantities.χZArray[idx_ZArray] *
                Weight_Matrix[idx_ZArray, idx_ZArrayInt]
            end
        end
    end
    LensingFunction.LensingEfficiencyArray .*=
    (last(CosmologicalGrid.ZArray)-first(CosmologicalGrid.ZArray))/
    (length(CosmologicalGrid.ZArray)-1)
end

"""
    ComputeWeightFunction(z::Float64, i::Int64,ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity, InstrumentResponse::InstrumentResponse,
    Cosmology::AbstractCosmology, CosmologicalGrid::CosmologicalGrid,
    LensingFunction::WLWeightFunction)

This function returns the Weak Lensing Weight Function, for a given redshift
``z`` and tomographic bin ``i``.
"""
function ComputeWeightFunction(z::Float64, i::Int64,
    ConvolvedDensity::AbstractConvolvedDensity, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, Cosmology::AbstractCosmology,
    CosmologicalGrid::CosmologicalGrid, LensingFunction::WLWeightFunction)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    return 1.5 * ComputeLensingEfficiency(z, i, ConvolvedDensity,
    AnalitycalDensity,InstrumentResponse, Cosmology, CosmologicalGrid,
    LensingFunction) * (Cosmology.H0/c_0)^2 * Cosmology.ΩM *
    (1. + z) * Computeχ(z, Cosmology)^2
end


"""
    ComputeWeightFunctionGrid!(LensingFunction::WLWeightFunction,
    ConvolvedDensity::AbstractConvolvedDensity, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities, Cosmology::AbstractCosmology)

This function evaluates the Weak Lensing Weight Function over the ``z``
grid and for all tomographic bins ``i``.
"""
function ComputeWeightFunctionGrid!(
    LensingFunction::WLWeightFunction,
    ConvolvedDensity::AbstractConvolvedDensity,
    CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    Cosmology::AbstractCosmology)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            LensingFunction.WeightFunctionArray[idx_ZBinArray, idx_ZArray] =
            1.5 * (Cosmology.H0/c_0)^2 * Cosmology.ΩM *
            (1. + CosmologicalGrid.ZArray[idx_ZArray]) *
            BackgroundQuantities.χZArray[idx_ZArray]^2 *
            LensingFunction.LensingEfficiencyArray[idx_ZBinArray, idx_ZArray] +
            LensingFunction.IntrinsicAlignmentArray[idx_ZBinArray, idx_ZArray]
        end
    end
end
