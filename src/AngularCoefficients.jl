abstract type IntegrationMethod end
struct NumericalIntegrationSimpson <: IntegrationMethod end
struct CustomTrapz <: IntegrationMethod end
struct BeyondLimber <: IntegrationMethod end

"""
    ComputeCℓ!(AngularCoefficients::AngularCoefficients,
    WeightFunctionA::AbstractWeightFunction,
    WeightFunctionB::AbstractWeightFunction,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::AbstractCosmology, CosmologicalGrid::CosmologicalGrid,
    PowerSpectrum::PowerSpectrum, ::NumericalIntegrationSimpson)

This function evaluates the Angular Coefficients for all tomographic bins and
multipole values. In order to evaluate the numerical integrals, it is used the
Simpson numerical method from
[NumericalIntegration.jl](https://github.com/dextorious/NumericalIntegration.jl)
"""
function  ComputeCℓ!(Cℓ::AbstractCℓ,
    WeightFunctionA::AbstractWeightFunction,
    WeightFunctionB::AbstractWeightFunction,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::AbstractCosmology, CosmologicalGrid::AbstractCosmologicalGrid,
    PowerSpectrum::AbstractPowerSpectrum, ::NumericalIntegrationSimpson)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for idx_a in 1:size(Cℓ.CℓArray, 2)
        for idx_b in idx_a:size(Cℓ.CℓArray, 3)
            for idx_l in 1:size(Cℓ.CℓArray, 1)
                Integrand = c_0 .*
                WeightFunctionA.WeightFunctionArray[idx_a, :] .*
                WeightFunctionB.WeightFunctionArray[idx_b, :] ./
                (BackgroundQuantities.HZArray .*
                BackgroundQuantities.rZArray.^2) .*
                PowerSpectrum.InterpolatedPowerSpectrum[idx_l,:]
                Cℓ.CℓArray[idx_l, idx_a,
                idx_b] = NumericalIntegration.integrate(CosmologicalGrid.ZArray,
                Integrand, SimpsonEvenFast())
                Cℓ.CℓArray[idx_l, idx_b,
                idx_a] = Cℓ.CℓArray[idx_l,
                idx_a, idx_b]
            end
        end
    end
end

"""
    ComputeCℓ!(AngularCoefficients::AngularCoefficients,
    WeightFunctionA::AbstractWeightFunction,
    WeightFunctionB::AbstractWeightFunction,
    BackgroundQuantities::BackgroundQuantities,
    w0waCDMCosmology::AbstractCosmology, CosmologicalGrid::CosmologicalGrid,
    PowerSpectrum::PowerSpectrum, ::CustomTrapz)

This function evaluates the Angular Coefficients for all tomographic bins and
multipole values. In order to evaluate the numerical integrals, it has been
implemented the Simpson rule.
"""
function  ComputeCℓ!(AngularCoefficients::AbstractCℓ,
    WeightFunctionA::AbstractWeightFunction,
    WeightFunctionB::AbstractWeightFunction,
    BackgroundQuantities::AbstractBackgroundQuantities,
    w0waCDMCosmology::AbstractCosmology, CosmologicalGrid::AbstractCosmologicalGrid,
    PowerSpectrum::AbstractPowerSpectrum, ::CustomTrapz)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    check = true
    while check == true
        Integrand = similar(AngularCoefficients.CℓArray) .*0
        Simpson_weights = SimpsonWeightArray(length(CosmologicalGrid.ZArray))
        @avx for i ∈ axes(AngularCoefficients.CℓArray,2),
            j ∈ axes(AngularCoefficients.CℓArray,3),
            l ∈ axes(AngularCoefficients.CℓArray,1)
            for z ∈ axes(CosmologicalGrid.ZArray,1)
                Integrand[l,i,j] += c_0 *
                WeightFunctionA.WeightFunctionArray[i, z] *
                WeightFunctionB.WeightFunctionArray[j, z] /
                (BackgroundQuantities.HZArray[z] *
                BackgroundQuantities.rZArray[z]^2) *
                PowerSpectrum.InterpolatedPowerSpectrum[l,z] *
                Simpson_weights[z]
            end
        end
        Integrand .*= (last(CosmologicalGrid.ZArray)-
        first(CosmologicalGrid.ZArray))/(length(CosmologicalGrid.ZArray)-1)
        AngularCoefficients.CℓArray = Integrand
        if any(isnan,Integrand)

        else
            check = false
        end
    end
end

function  ComputeCℓ!(Cℓ::AbstractCℓ,
    TransferFunctionA::AbstractTransferFunction,
    TransferFunctionB::AbstractTransferFunction,
    w0waCDMCosmology::AbstractCosmology, CosmologicalGrid::AbstractCosmologicalGrid,
    PowerSpectrum::AbstractPowerSpectrum, ::BeyondLimber)
    Integrand = zeros(size(Cℓ.CℓArray))
    WeightsMatrix =
    UnevenTrapzWeightMatrix(CosmologicalGrid.KBeyondLimberArray)
    @avx for i ∈ axes(Cℓ.CℓArray,2),
        j ∈ axes(Cℓ.CℓArray,3),
        l ∈ axes(Cℓ.CℓArray,1)
        for k ∈ axes(CosmologicalGrid.KBeyondLimberArray,2)
            Integrand[l,i,j] += 2 / π * WeightsMatrix[l, k] *
            TransferFunctionA.TransferFunctionArray[i, l, k] *
            TransferFunctionB.TransferFunctionArray[j, l, k] *
            PowerSpectrum.InterpolatedPowerSpectrumBeyondLimber[l,k] *
            CosmologicalGrid.KBeyondLimberArray[l, k]^2
        end
    end
    Cℓ.CℓArray = Integrand
end
