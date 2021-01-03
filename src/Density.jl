abstract type AsbtractDensity end
abstract type AnalitycalDensity <: AsbtractDensity end
abstract type ConvolvedDensity <: AsbtractDensity end
abstract type InstrumentResponse end




"""
    AnalitycalDensityStruct(z0::Float64 = 0.9/sqrt(2.), zmin::Float64 = 0.001,
    zmax::Float64 = 2.5, surfacedensity::Float64 = 30.,
    normalization::Float64 = 1.)

This struct contains the parameters of the source galaxy density as given by the [official Euclid forecast](https://arxiv.org/abs/1910.09273), whose expression is given by:
```math
n(z)\\propto\\left(\\frac{z}{z_0}\\right)^2 \\exp{\\left(-\\left(\\frac{z}{z_0}\\right)^{-3/2}\\right)}
```
The parameters contained in this struct are
- ``z_{min}`` and ``z_{max}``, the minimum and the maximum redshift considered

- ``z_0``, the parameter present in the galaxy distribution

- surfacedensity , the value of the galaxy source density integrated between ``z_{min}`` and ``z_{max}``

- normalization, the value of parameter which multiplies the source dennsity in order to match the correct surface density
"""
@kwdef mutable struct AnalitycalDensityStruct <: AnalitycalDensity
    z0::Float64 = 0.9/sqrt(2.)
    zmin::Float64 = 0.001
    zmax::Float64 = 2.5
    surfacedensity::Float64 = 30.
    normalization::Float64 = 1.
end

"""
    ConvolvedDensityStruct(AnalitycalDensity::AnalitycalDensity = AnalitycalDensityStruct(),
    InstrumentResponse::InstrumentResponse = InstrumentResponseStruct()
    zbinarray::Vector{Float64} = Array([0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.50])
    densityarraynormalization::Vector{Float64} = ones(length(zbinarray)-1)
    densitygridarray::AbstractArray{Float64, 2} = ones(length(zbinarray)-1, 300))

This function evaluate the normalization of the convolved density function.
"""
@kwdef mutable struct ConvolvedDensityStruct <: ConvolvedDensity
    AnalitycalDensity::AnalitycalDensity = AnalitycalDensityStruct()
    InstrumentResponse::InstrumentResponse = InstrumentResponseStruct()
    zbinarray::Vector{Float64} = Array([0.001, 0.418, 0.560, 0.678, 0.789,
    0.900, 1.019, 1.155, 1.324, 1.576, 2.50])
    densityarraynormalization::Vector{Float64} = ones(length(zbinarray)-1)
    densitygridarray::AbstractArray{Float64, 2} = ones(length(zbinarray)-1, 300)
end

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

This function normalize AnalitycalDensityStruct in order to have the correct value of the surface density once integrated.
"""
function NormalizeAnalitycalDensityStruct(densityparameters::AnalitycalDensity)
    int, err = QuadGK.quadgk(x -> ComputeDensityFunction(x, densityparameters),
    densityparameters.zmin, densityparameters.zmax, rtol=1e-12)
    densityparameters.normalization *= (densityparameters.surfacedensity/int)
end

"""
    InstrumentResponseStruct(cb::Float64 = 1.0, zb::Float64 = 0.0,
    σb::Float64 = 0.05, co::Float64 = 1.0, zo::Float64 = 0.1,
    σo::Float64 = 0.05, fout::Float64 = 0.1)

When we measure the redshift of a galaxy with redshit ``z``, we will measure a
redshift ``z_p`` with a probability given by
[the following expression](https://arxiv.org/abs/1910.09273):
```math
p(z_p|z)  = \\frac{1-f_{out}}{\\sqrt{2 \\pi} \\sigma_{b}(1+z)} \\exp \\left(
-\\frac{1}{2}\\left(\\frac{z-c_{b} z_{b}-z_{b}}{\\sigma_{b}(1+z)}\\right)^{2}
\\right) + \\frac{f_{out}}{\\sqrt{2 \\pi} \\sigma_{\\mathrm{o}}(1+z)} \\exp
\\left(-\\frac{1}{2}\\left(\\frac{z-c_{o} z_{p}-z_{o}}{\\sigma_{o}(1+z)}
\\right)^{2}\\right)
```
This struct contains all these parameters.
"""
@kwdef struct InstrumentResponseStruct <: InstrumentResponse
    cb::Float64   = 1.0
    zb::Float64   = 0.0
    σb::Float64   = 0.05
    co::Float64   = 1.0
    zo::Float64   = 0.1
    σo::Float64   = 0.05
    fout::Float64 = 0.1
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
    ComputeDensityFunction(z::Float64, i::Int64, convolveddensity::ConvolvedDensity)

This function computes the Convolved density function for a single bin at a given redshift ``z``.
"""
function ComputeConvolvedDensityFunction(z::Float64, i::Int64,
    convolveddensity::ConvolvedDensity)
    int, err = QuadGK.quadgk(x -> ComputeInstrumentResponse(z, x,
    convolveddensity.InstrumentResponse),
    convolveddensity.zbinarray[i],
    convolveddensity.zbinarray[i+1], rtol=1e-12)
    return int*ComputeDensityFunction(z, convolveddensity.AnalitycalDensity)*
    convolveddensity.densityarraynormalization[i]
end

"""
    NormalizeConvolvedDensityStruct(convolveddensity::ConvolvedDensity)

This function normalizes ConvolvedDensity such that the integrals of the convolved densities are normalized to 1.
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
    ComputeDensityFunction(CosmoGrid::CosmoGrid, ConvolvedDensityStruct::ConvolvedDensity)

This function computes the convolved density function for all tomographic bins
on the ``z``-grid provided by CosmoGrid.
"""
function ComputeDensityFunctionConvolvedGrid(CosmoGrid::CosmoGrid,
    ConvolvedDensityStruct::ConvolvedDensity)
    ConvolvedDensityStruct.densitygridarray = zeros(Float64,length(ConvolvedDensityStruct.zbinarray)-1,
    length(CosmoGrid.zgrid))
    for idx_zbinarray in 1:length(ConvolvedDensityStruct.zbinarray)-1
        for idx_zgrid in 1:length(CosmoGrid.zgrid)
            ConvolvedDensityStruct.densitygridarray[idx_zbinarray, idx_zgrid] =
            ComputeConvolvedDensityFunction(CosmoGrid.zgrid[idx_zgrid],
            idx_zbinarray, ConvolvedDensityStruct)
        end
    end
end
