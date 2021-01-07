"""
    ComputeDensityFunction(z::Float64, AnalitycalDensity::AnalitycalDensity = AnalitycalDensityStruct())

This function returns the source density for a given redshift ``z``.
"""
function ComputeDensityFunction(z::Float64, densityparameters::AnalitycalDensity)
    return (z/densityparameters.z0)^2*exp(-(z/densityparameters.z0)^(3. / 2.))*
    densityparameters.normalization
end

"""
    NormalizeAnalitycalDensityStruct(densityparameters::AnalitycalDensity)

This function normalize AnalitycalDensityStruct in order to have the correct
value of the surface density once integrated.
"""
function NormalizeAnalitycalDensityStruct(densityparameters::AnalitycalDensity)
    int, err = QuadGK.quadgk(x -> ComputeDensityFunction(x, densityparameters),
    densityparameters.zmin, densityparameters.zmax, rtol=1e-12)
    densityparameters.normalization *= (densityparameters.surfacedensity/int)
end


"""
    ComputeInstrumentResponse(z::Float64, zp::Float64, instrumentresponse::InstrumentResponse)

This function computes the probability that we actually measure a redshift
``z_p`` if the real redshift is ``z``.
"""
function ComputeInstrumentResponse(z::Float64, zp::Float64,
    instrumentresponse::InstrumentResponseStruct)
    prob_z = (1-instrumentresponse.fout)/(sqrt(2*pi)*instrumentresponse.σb*(1+z))*
    exp(-0.5*((z-instrumentresponse.cb*zp-instrumentresponse.zb)/(instrumentresponse.σb*(1+z)))^2)
    +instrumentresponse.fout/(sqrt(2*pi)*instrumentresponse.σo*(1+z))*
    exp(-0.5*((z-instrumentresponse.co*zp-instrumentresponse.zo)/(instrumentresponse.σo*(1+z)))^2)
    return prob_z
end

"""
    ComputeConvolvedDensityFunction(z::Float64, i::Int64, convolveddensity::ConvolvedDensity)

This function computes the Convolved density function for a single bin at a
given redshift ``z``.
"""
function ComputeConvolvedDensityFunction(z::Float64, i::Int64,
    convolveddensity::ConvolvedDensity)
    int, err = QuadGK.quadgk(x -> ComputeInstrumentResponse(z, x,
    convolveddensity.InstrumentResponse),
    convolveddensity.ZBinArray[i],
    convolveddensity.ZBinArray[i+1], rtol=1e-12)
    return int*ComputeDensityFunction(z, convolveddensity.AnalitycalDensity)*
    convolveddensity.densityarraynormalization[i]
end

"""
    NormalizeConvolvedDensityStruct(convolveddensity::ConvolvedDensity)

This function normalizes ConvolvedDensity such that the integrals of the
convolved densities are normalized to 1.
"""
function NormalizeConvolvedDensityStruct(convolveddensity::ConvolvedDensity)
    for idx in 1:length(convolveddensity.densityarraynormalization)
        int, err = QuadGK.quadgk(x -> ComputeConvolvedDensityFunction(x, idx,
        convolveddensity),
        convolveddensity.AnalitycalDensity.zmin,
        convolveddensity.AnalitycalDensity.zmax, rtol=1e-12)
        convolveddensity.densityarraynormalization[idx] /= int
    end
end

"""
    ComputeConvolvedDensityFunctionGrid(CosmologicalGrid::CosmologicalGrid, ConvolvedDensityStruct::ConvolvedDensity)

This function computes the convolved density function for all tomographic bins
on the ``z``-grid provided by CosmologicalGrid.
"""
function ComputeConvolvedDensityFunctionGrid(CosmologicalGrid::CosmologicalGrid,
    ConvolvedDensityStruct::ConvolvedDensity)
    ConvolvedDensityStruct.densitygridarray = zeros(Float64,length(ConvolvedDensityStruct.ZBinArray)-1,
    length(CosmologicalGrid.ZArray))
    for idx_ZBinArray in 1:length(ConvolvedDensityStruct.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            ConvolvedDensityStruct.densitygridarray[idx_ZBinArray, idx_ZArray] =
            ComputeConvolvedDensityFunction(CosmologicalGrid.ZArray[idx_ZArray],
            idx_ZBinArray, ConvolvedDensityStruct)
        end
    end
end
