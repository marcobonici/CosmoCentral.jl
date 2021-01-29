abstract type IntegrationMethod end

struct NumericalIntegrationSimpson <: IntegrationMethod end
struct CustomTrapz <: IntegrationMethod end


"""
    ComputeAngularCoefficients(AngularCoefficients::AngularCoefficients,
    WeightFunctionA::AbstractWeightFunction,
    WeightFunctionB::AbstractWeightFunction,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::AbstractCosmology, CosmologicalGrid::CosmologicalGrid,
    PowerSpectrum::PowerSpectrum, ::IntegrationMethod)

This function evaluates the Angular Coefficients for all tomographic bins and
multipole values.
"""
function  ComputeAngularCoefficients(AngularCoefficients::AngularCoefficients,
    WeightFunctionA::AbstractWeightFunction,
    WeightFunctionB::AbstractWeightFunction,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::AbstractCosmology, CosmologicalGrid::CosmologicalGrid,
    PowerSpectrum::PowerSpectrum, ::NumericalIntegrationSimpson)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for idx_a in 1:size(AngularCoefficients.AngularCoefficientsArray, 2)
        for idx_b in idx_a:size(AngularCoefficients.AngularCoefficientsArray, 3)
            for idx_l in 1:size(AngularCoefficients.AngularCoefficientsArray, 1)
                Integrand = c_0 .*
                WeightFunctionA.WeightFunctionArray[idx_a, :] .*
                WeightFunctionB.WeightFunctionArray[idx_b, :] ./
                (BackgroundQuantities.HZArray .*
                BackgroundQuantities.rZArray.^2) .*
                PowerSpectrum.InterpolatedPowerSpectrum[idx_l,:]
                AngularCoefficients.AngularCoefficientsArray[idx_l, idx_a,
                idx_b] = NumericalIntegration.integrate(CosmologicalGrid.ZArray,
                Integrand, SimpsonEvenFast())
                AngularCoefficients.AngularCoefficientsArray[idx_l, idx_b,
                idx_a] = AngularCoefficients.AngularCoefficientsArray[idx_l,
                idx_a, idx_b]
            end
        end
    end
end

function  ComputeAngularCoefficients(AngularCoefficients::AngularCoefficients,
    WeightFunctionA::AbstractWeightFunction,
    WeightFunctionB::AbstractWeightFunction,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::AbstractCosmology, CosmologicalGrid::CosmologicalGrid,
    PowerSpectrum::PowerSpectrum, ::CustomTrapz)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    Integrand = zeros(size(AngularCoefficients.AngularCoefficientsArray,1),
    size(AngularCoefficients.AngularCoefficientsArray,2),
    size(AngularCoefficients.AngularCoefficientsArray,3))
    @avx for i ∈ axes(AngularCoefficients.AngularCoefficientsArray,2),
        j ∈ axes(AngularCoefficients.AngularCoefficientsArray,3),
        l ∈ axes(AngularCoefficients.AngularCoefficientsArray,1)
        for z ∈ axes(CosmologicalGrid.ZArray,1)
            Integrand[l,i,j] += c_0 *
            WeightFunctionA.WeightFunctionArray[i, z] *
            WeightFunctionB.WeightFunctionArray[j, z] /
            (BackgroundQuantities.HZArray[z] *
            BackgroundQuantities.rZArray[z]^2) *
            PowerSpectrum.InterpolatedPowerSpectrum[l,z]
        end
    end
    Integrand .*= (last(CosmologicalGrid.ZArray)-first(CosmologicalGrid.ZArray))/(length(CosmologicalGrid.ZArray)-1)
    AngularCoefficients.AngularCoefficientsArray = Integrand
end
