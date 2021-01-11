"""
    ComputeDensityFunction(z::Float64, AnalitycalDensity::AnalitycalDensity = AnalitycalDensityStruct())

This function returns the source density for a given redshift ``z``.
"""
function ComputeDensityFunction(z::Float64, AnalitycalDensity::AnalitycalDensity)
    return ((z/AnalitycalDensity.Z0)^2)*exp(-(z/AnalitycalDensity.Z0)^(3. / 2.))*
    AnalitycalDensity.Normalization
end

"""
    NormalizeAnalitycalDensityStruct(densityparameters::AnalitycalDensity)

This function normalize AnalitycalDensityStruct in order to have the correct
value of the surface density once integrated.
"""
function NormalizeAnalitycalDensityStruct(AnalitycalDensity::AnalitycalDensity)
    int, err = QuadGK.quadgk(x -> ComputeDensityFunction(x, AnalitycalDensity),
    AnalitycalDensity.ZMin, AnalitycalDensity.ZMax, rtol=1e-12)
    AnalitycalDensity.Normalization *= (AnalitycalDensity.SurfaceDensity/int)
end


"""
    ComputeInstrumentResponse(z::Float64, zp::Float64, instrumentresponse::InstrumentResponse)

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
    ComputeConvolvedDensityFunction(z::Float64, i::Int64, convolveddensity::ConvolvedDensity)

This function computes the Convolved density function for a single bin at a
given redshift ``z``.
"""
function ComputeConvolvedDensityFunction(z::Float64, i::Int64,
    ConvolvedDensity::ConvolvedDensity, AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse)
    int, err = QuadGK.quadgk(x -> ComputeInstrumentResponse(z, x,
    InstrumentResponse),
    ConvolvedDensity.ZBinArray[i],
    ConvolvedDensity.ZBinArray[i+1], rtol=1e-12)
    return int*ComputeDensityFunction(z, AnalitycalDensity)*
    ConvolvedDensity.DensityNormalizationArray[i]
end

"""
    NormalizeConvolvedDensityStruct(convolveddensity::ConvolvedDensity)

This function normalizes ConvolvedDensity such that the integrals of the
convolved densities are normalized to 1.
"""
function NormalizeConvolvedDensityStruct(ConvolvedDensity::ConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse,
    CosmologicalGrid::CosmologicalGrid)
    for idx in 1:length(ConvolvedDensity.DensityNormalizationArray)
        int, err = QuadGK.quadgk(x -> ComputeConvolvedDensityFunction(x, idx,
        ConvolvedDensity, AnalitycalDensity, InstrumentResponse),
        first(CosmologicalGrid.ZArray),
        last(CosmologicalGrid.ZArray), rtol=1e-12)
        ConvolvedDensity.DensityNormalizationArray[idx] /= int
    end
end

"""
    ComputeConvolvedDensityFunctionGrid(CosmologicalGrid::CosmologicalGrid, ConvolvedDensityStruct::ConvolvedDensity)

This function computes the convolved density function for all tomographic bins
on the ``z``-grid provided by CosmologicalGrid.
"""
function ComputeConvolvedDensityFunctionGrid(CosmologicalGrid::CosmologicalGrid,
    ConvolvedDensityStruct::ConvolvedDensity,
    AnalitycalDensity::AnalitycalDensity,
    InstrumentResponse::InstrumentResponse)
    ConvolvedDensityStruct.DensityGridArray = zeros(Float64,length(ConvolvedDensityStruct.ZBinArray)-1,
    length(CosmologicalGrid.ZArray))
    for idx_ZBinArray in 1:length(ConvolvedDensityStruct.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            ConvolvedDensityStruct.DensityGridArray[idx_ZBinArray, idx_ZArray] =
            ComputeConvolvedDensityFunction(CosmologicalGrid.ZArray[idx_ZArray],
            idx_ZBinArray, ConvolvedDensityStruct, AnalitycalDensity,
            InstrumentResponse)
        end
    end
end
