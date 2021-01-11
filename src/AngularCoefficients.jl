"""
    ComputeAngularCoefficients(AngularCoefficients::AngularCoefficients,
    WeightFunctionA::GCWeightFunction, WeightFunctionB::GCWeightFunction,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology, CosmologicalGrid::CosmologicalGrid,
    PowerSpectrum::PowerSpectrum)

This function evaluates the Angular Coefficients for all tomographic bins and
multipole values.
"""
function  ComputeAngularCoefficients(AngularCoefficients::AngularCoefficients,
    WeightFunctionA::GCWeightFunction, WeightFunctionB::GCWeightFunction,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::w0waCDMCosmology, CosmologicalGrid::CosmologicalGrid,
    PowerSpectrum::PowerSpectrum)
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
            end
        end
    end
end
