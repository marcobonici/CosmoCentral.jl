"""
    ComputeDensityFunction(z::Float64, AnalitycalDensity::AnalitycalDensity)

This function returns the source density for a given redshift ``z``.
"""
function ComputeDensity(z::T, AnalitycalDensity::AnalitycalDensity) where T
    return ((z/AnalitycalDensity.Z0)^2)*exp(-(z/AnalitycalDensity.Z0)^(3. / 2.))*
    AnalitycalDensity.Normalization
end

"""
    NormalizeAnalitycalDensity!(AnalitycalDensity::AnalitycalDensity)

This function normalize AnalitycalDensity in order to have the correct
value of the surface density once integrated.
"""
function NormalizeAnalitycalDensity!(AnalitycalDensity::AnalitycalDensity)
    int, err = quadgk(x -> ComputeDensity(x, AnalitycalDensity),
    AnalitycalDensity.ZMin, AnalitycalDensity.ZMax, rtol=1e-12)
    AnalitycalDensity.Normalization *= (AnalitycalDensity.SurfaceDensity/int)
end


"""
    ComputeInstrumentResponse(z::Float64, zp::Float64,
    InstrumentResponse::InstrumentResponse)

This function computes the probability that we actually measure a redshift
``z_p`` if the real redshift is ``z``.
"""
function ComputeInstrumentResponse(z::T, zp::T,
    InstrumentResponse::InstrumentResponse) where T
    prob_z =
    (1-InstrumentResponse.fout)/(sqrt(2*pi)*InstrumentResponse.σb*(1+z))*
    exp(-0.5*((z-InstrumentResponse.cb*zp-InstrumentResponse.zb)/
    (InstrumentResponse.σb*(1+z)))^2)+
    InstrumentResponse.fout/(sqrt(2*pi)*InstrumentResponse.σo*(1+z))*
    exp(-0.5*((z-InstrumentResponse.co*zp-InstrumentResponse.zo)/
    (InstrumentResponse.σo*(1+z)))^2)
    return prob_z
end

"""
    ComputeConvolvedDensity(z::Float64, i::Int64, ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity, InstrumentResponse::InstrumentResponse)

This function computes the Convolved density function for a single bin at a
given redshift ``z``.
"""
function ComputeConvolvedDensity(z::T, i::I,
    ConvolvedDensity::AbstractConvolvedDensity, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse) where {T, I}
    int, err = quadgk(x -> ComputeInstrumentResponse(z, x,
    InstrumentResponse),
    ConvolvedDensity.ZBinArray[i],
    ConvolvedDensity.ZBinArray[i+1], rtol=1e-12)
    return int*ComputeDensity(z,
    AnalitycalDensity) * ConvolvedDensity.DensityNormalizationArray[i]
end

"""
    NormalizeConvolvedDensity!(ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity,InstrumentResponse::InstrumentResponse,
    CosmologicalGrid::CosmologicalGrid)

This function normalizes ConvolvedDensity such that the integrals of the
convolved densities are normalized to 1.
"""
function NormalizeConvolvedDensity!(ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity, InstrumentResponse::InstrumentResponse,
    CosmologicalGrid::CosmologicalGrid)
    for idx in 1:length(ConvolvedDensity.DensityNormalizationArray)
        int, err = quadgk(x -> ComputeConvolvedDensity(x, idx,
        ConvolvedDensity, AnalitycalDensity, InstrumentResponse),
        first(CosmologicalGrid.ZArray),
        last(CosmologicalGrid.ZArray), rtol=1e-12)
        ConvolvedDensity.DensityNormalizationArray[idx] /= int
    end
end

function ComputeSurfaceDensityBins!(ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity)
    for idx in 1:length(ConvolvedDensity.SurfaceDensityArray)
        int, err = quadgk(x -> ComputeDensity(x, AnalitycalDensity),
        ConvolvedDensity.ZBinArray[idx], ConvolvedDensity.ZBinArray[idx+1], rtol=1e-12)
        ConvolvedDensity.SurfaceDensityArray[idx] = int
    end
end


"""
    ComputeConvolvedDensityGrid!(CosmologicalGrid::CosmologicalGrid,
    ConvolvedDensity::AbstractConvolvedDensity, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse)

This function computes the convolved density function for all tomographic bins
on the ``z``-grid provided by CosmologicalGrid.
"""
function ComputeConvolvedDensityGrid!(CosmologicalGrid::CosmologicalGrid,
    ConvolvedDensity::AbstractConvolvedDensity, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse)
    ConvolvedDensity.DensityGridArray = zeros(Float64,length(ConvolvedDensity.ZBinArray)-1,
    length(CosmologicalGrid.ZArray))
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            ConvolvedDensity.DensityGridArray[idx_ZBinArray, idx_ZArray] =
            ComputeConvolvedDensity(CosmologicalGrid.ZArray[idx_ZArray],
            idx_ZBinArray, ConvolvedDensity, AnalitycalDensity,
            InstrumentResponse)
        end
    end
end

function ShiftConvolvedDensityGrid!(CosmologicalGrid::CosmologicalGrid,
    ConvolvedDensity::AbstractConvolvedDensity)
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        InterpConvDens = Dierckx.Spline1D(CosmologicalGrid.ZArray,
        ConvolvedDensity.DensityGridArray[idx_ZBinArray,:], k=4,
        bc="extrapolate")
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            ConvolvedDensity.DensityGridArray[idx_ZBinArray,idx_ZArray] =
            InterpConvDens(CosmologicalGrid.ZArray[idx_ZArray]+
            ConvolvedDensity.ShiftArray[idx_ZBinArray])
        end
    end
end

function CreateAnalitycalDensity(DensDict::Dict)
    andensity = AnalitycalDensity()
    andensity.Z0 = DensDict["Z0"]
    andensity.ZMin = DensDict["ZMin"]
    andensity.ZMax = DensDict["ZMax"]
    andensity.SurfaceDensity = DensDict["SurfaceDensity"]
    return andensity
end

function CreateDensity(DensityDict::Dict, CosmologicalGrid::CosmologicalGrid)
    if DensityDict["model"] == "AnalitycalDensity"
        analyticaldensity = CreateAnalitycalDensity(DensityDict["parameters"])
        instrumentresponse = InstrumentResponse()
        convolveddensity = ConvolvedDensity(DensityGridArray =
        ones(10, length(CosmologicalGrid.ZArray)))
        NormalizeConvolvedDensity!(convolveddensity, analyticaldensity,
        instrumentresponse, CosmologicalGrid)
        ComputeConvolvedDensityGrid!(CosmologicalGrid, convolveddensity,
        analyticaldensity, instrumentresponse)
    elseif DensityDict["model"] == "KinematicLensing"
        nigz = npzread("/home/mbonici/Downloads/Materiale/spectroscopic/13bin/n_i.npy")
        bar_niz = npzread("/home/mbonici/Downloads/Materiale/spectroscopic/13bin/bar_n_i_KL.npy") .* 8.461595013198466e-8#horrible patch
        ZArray = npzread("/home/mbonici/Downloads/Materiale/spectroscopic/13bin/z.npy")
        zbinarray = [0.901, 0.97, 1.039, 1.108, 1.177, 1.247, 1.316, 1.385, 1.454, 1.524, 1.593, 1.662, 1.731, 1.799]
        convolveddensity = ConvolvedDensity(DensityGridArray = nigz,
        DensityNormalizationArray = ones(13), ZBinArray = zbinarray,
        SurfaceDensityArray = bar_niz)
    end
    return convolveddensity
end
