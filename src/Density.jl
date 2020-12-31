abstract type AsbtractDensity end
abstract type AnalitycalDensity <: AsbtractDensity end

@kwdef mutable struct AnalitycalDensityStruct <: AnalitycalDensity
    z0::Float64 = 0.9/sqrt(2.)
    zmin::Float64 = 0.001
    zmax::Float64 = 2.5
    surfacedensity::Float64 = 30
    normalization::Float64 = 1
end

function ComputeDensityFunction(z::Float64, densityparameters::AnalitycalDensity)
    return (z/densityparameters.z0)^2*exp(-(z/densityparameters.z0)^(3. / 2.))/
    densityparameters.normalization
end

function NormalizeAnalitycalDensityStruct(densityparameters::AnalitycalDensity)
    densityparameters.normalization = 1
    normalization, err = QuadGK.quadgk(x -> ComputeDensityFunction(x, densityparameters),
    densityparameters.zmin, densityparameters.zmax, rtol=1e-12)
    densityparameters.normalization = normalization*densityparameters.surfacedensity
    return densityparameters
end
