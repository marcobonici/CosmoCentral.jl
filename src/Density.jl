abstract type AsbtractDensity end
abstract type AnalitycalDensity <: AsbtractDensity end
abstract type ConvolvedDensity end
abstract type InstrumentResponse end

@kwdef mutable struct AnalitycalDensityStruct <: AnalitycalDensity
    z0::Float64 = 0.9/sqrt(2.)
    zmin::Float64 = 0.001
    zmax::Float64 = 2.5
    surfacedensity::Float64 = 30.
    normalization::Float64 = 1.
end

"""
ComputeDensityFunction(``z``, params)

This function returns the source density for a given redshift ``z``
"""
function ComputeDensityFunction(z::Float64, densityparameters::AnalitycalDensity)
    return (z/densityparameters.z0)^2*exp(-(z/densityparameters.z0)^(3. / 2.))*
    densityparameters.normalization
end

"""
ComputeDensityFunction(params)

This function modifies the normalization constant in the AnalitycalDensityStruct in order to have the same value of the surface density once integrated.
"""
function NormalizeAnalitycalDensityStruct(densityparameters::AnalitycalDensity)
    int, err = QuadGK.quadgk(x -> ComputeDensityFunction(x, densityparameters),
    densityparameters.zmin, densityparameters.zmax, rtol=1e-12)
    densityparameters.normalization *= (densityparameters.surfacedensity/int)
    return densityparameters
end

@kwdef struct InstrumentResponseStruct <: InstrumentResponse
    cb::Float64   = 1.0
    zb::Float64   = 0.0
    σb::Float64   = 0.05
    co::Float64   = 1.0
    zo::Float64   = 0.1
    σ0::Float64   = 0.05
    fout::Float64 = 0.1
end


"""
ComputeInstrumentResponse(``z``, params)

This function returns the instrument response.
"""
function ComputeInstrumentResponse(z::Float64, zp::Float64,
    instrumentresponse::InstrumentResponseStruct)
    prob_z = (1-instrumentresponse.fout)/(sqrt(2*pi)*instrumentresponse.σb*(1+z))*
    exp(-0.5*((z-instrumentresponse.cb*zp-instrumentresponse.zb)/(instrumentresponse.σb*(1+z)))^2)
    +instrumentresponse.fout/(sqrt(2*pi)*instrumentresponse.σo*(1+z))*
    exp(-0.5*((z-instrumentresponse.co*zp-instrumentresponse.zo)/(instrumentresponse.σo*(1+z)))^2)
    return prob_z
end
