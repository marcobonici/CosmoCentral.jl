"""
    ComputeDensityFunction(z::Float64,
    AnalitycalDensity::AnalitycalDensityStruct)

This function returns the source density for a given redshift ``z``.
"""
function ComputeDensity(z::Float64,
    AnalitycalDensity::AnalitycalDensity)
    return ((z/AnalitycalDensity.Z0)^2)*exp(-(z/AnalitycalDensity.Z0)^(3. / 2.))*
    AnalitycalDensity.Normalization
end

"""
    NormalizeAnalitycalDensity!(densityparameters::AnalitycalDensityStruct)

This function normalize AnalitycalDensityStruct in order to have the correct
value of the surface density once integrated.
"""
function NormalizeAnalitycalDensity!(
    AnalitycalDensity::AnalitycalDensity)
    int, err = QuadGK.quadgk(x -> ComputeDensity(x, AnalitycalDensity),
    AnalitycalDensity.ZMin, AnalitycalDensity.ZMax, rtol=1e-12)
    AnalitycalDensity.Normalization *= (AnalitycalDensity.SurfaceDensity/int)
end


"""
    ComputeInstrumentResponse(z::Float64, zp::Float64,
    InstrumentResponse::InstrumentResponse)

This function computes the probability that we actually measure a redshift
``z_p`` if the real redshift is ``z``.
"""
function ComputeInstrumentResponse(z::Float64, zp::Float64,
    InstrumentResponse::InstrumentResponse)
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
    ComputeConvolvedDensity(z::Float64, i::Int64,
    ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse)

This function computes the Convolved density function for a single bin at a
given redshift ``z``.
"""
function ComputeConvolvedDensity(z::Float64, i::Int64,
    ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse)
    int, err = QuadGK.quadgk(x -> ComputeInstrumentResponse(z, x,
    InstrumentResponse),
    ConvolvedDensity.ZBinArray[i],
    ConvolvedDensity.ZBinArray[i+1], rtol=1e-12)
    return int*ComputeDensity(z,
    AnalitycalDensity) * ConvolvedDensity.DensityNormalizationArray[i]
end

"""
    NormalizeConvolvedDensity!(ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse, CosmologicalGrid::CosmologicalGrid)

This function normalizes ConvolvedDensity such that the integrals of the
convolved densities are normalized to 1.
"""
function NormalizeConvolvedDensity!(
    ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse, CosmologicalGrid::CosmologicalGrid)
    for idx in 1:length(ConvolvedDensity.DensityNormalizationArray)
        int, err = QuadGK.quadgk(x -> ComputeConvolvedDensity(x, idx,
        ConvolvedDensity, AnalitycalDensity, InstrumentResponse),
        first(CosmologicalGrid.ZArray),
        last(CosmologicalGrid.ZArray), rtol=1e-12)
        ConvolvedDensity.DensityNormalizationArray[idx] /= int
    end
end

"""
    ComputeConvolvedDensityGrid!(CosmologicalGrid::CosmologicalGrid,
    ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensityStruct,
    InstrumentResponse::InstrumentResponse)

This function computes the convolved density function for all tomographic bins
on the ``z``-grid provided by CosmologicalGrid.
"""
function ComputeConvolvedDensityGrid!(CosmologicalGrid::CosmologicalGrid,
    ConvolvedDensity::AbstractConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity,
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
