abstract type AsbtractDensity end
abstract type AnalitycalDensity <: AsbtractDensity end
abstract type ConvolvedDensity <: AsbtractDensity end
abstract type InstrumentResponse end

@kwdef mutable struct AnalitycalDensityStruct <: AnalitycalDensity
    z0::Float64 = 0.9/sqrt(2.)
    zmin::Float64 = 0.001
    zmax::Float64 = 2.5
    surfacedensity::Float64 = 30.
    normalization::Float64 = 1.
end

@kwdef mutable struct ConvolvedDensityStruct <: ConvolvedDensity
    AnalitycalDensity::AnalitycalDensity = AnalyticalDensityStruct()
    InstrumentResponse::InstrumentResponse = InstrumentResponseStruct()
    zbinarray::Vector{Float64} = Array([0.001, 0.418, 0.560, 0.678, 0.789,
    0.900, 1.019, 1.155, 1.324, 1.576, 2.50])
    densityarraynormalization::Vector{Float64} = ones(length(zbinarray)-1)
    densitygridarray::AbstractArray{Float64, 2} = ones(length(zbinarray)-1, 300)
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

function ComputeConvolvedDensityFunction(z::Float64, i::Int64,
    convolveddensity::ConvolvedDensity)
    int, err = QuadGK.quadgk(x -> ComputeInstrumentResponse(z, x,
    convolveddensity.InstrumentResponse),
    convolveddensity.analitycaldensity.zbinarray[i],
    convolveddensity.analitycaldensity.zbinarray[i+1], rtol=1e-12)
    return int*ComputeDensityFunction(z, convolveddensity.analitycaldensity)*
    convolveddensity.densityarraynormalization[i]
end

"""
ComputeDensityFunction(params)

This function modifies the normalization constant in the AnalitycalDensityStruct in order to have the same value of the surface density once integrated.
"""
function NormalizeConvolvedDensityStruct(convolveddensity::ConvolvedDensity)
    for idx in 1:length(convolveddensity.zarraynormalization)
        int, err = QuadGK.quadgk(x -> ComputeDensityFunction(x, convolveddensity),
        convolveddensity.zmin, convolveddensity.zmax, rtol=1e-12)
        convolveddensity.densityarraynormalization[i] /= int
    end
    return convolveddensity
end


function ComputeDensityFunctionConvolvedGrid(CosmoGrid::CosmoGrid,
    InstrumentResponse::InstrumentResponse, analitycaldensity::AnalitycalDensity)
    grid = zeros(Float64,length(densityparameters.zbinarray)-1, length(CosmoGrid.zgrid))
    for idx_zbinarray in 1:length(densityparameters.zbinarray)-1
        for idx_zgrid in 1:length(CosmoGrid.zgrid)
            grid[idx_zbinarray, idx_zgrid] = ComputeConvolvedDensityFunction(
            CosmoGrid.zgrid[idx_zgrid], analitycaldensity.zbinarray[idx_zbinarray],
            InstrumentResponse, densityparameters)
        end
        normalization = NumericalIntegration.integrate(CosmoGrid.zgrid,
        grid[idx_zbinarray,:], SimpsonEven())
        grid[idx_zbinarray,:] ./= normalization
    end
    return grid
end
