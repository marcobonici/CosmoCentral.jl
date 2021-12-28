abstract type IntegrationMethod end
struct NumericalIntegrationSimpson <: IntegrationMethod end
struct CustomSimpson <: IntegrationMethod end
struct BeyondLimber <: IntegrationMethod end

"""
    ComputeCℓ!(Cℓ::AbstractCℓ, WeightFunctionA::AbstractWeightFunction,
    WeightFunctionB::AbstractWeightFunction, BackgroundQuantities::BackgroundQuantities,
    ::AbstractCosmology, CosmologicalGrid::AbstractCosmologicalGrid,
    PowerSpectrum::AbstractPowerSpectrum, ::NumericalIntegrationSimpson)

This function evaluates the Angular Coefficients for all tomographic bins and
multipole values. In order to evaluate the numerical integrals, it is employed the
Simpson numerical method from
[NumericalIntegration.jl](https://github.com/dextorious/NumericalIntegration.jl) . This is
not the fastest method available, but can be used as a benchmark to check consistency.
"""
function  ComputeCℓ!(Cℓ::AbstractCℓ, WeightFunctionA::AbstractWeightFunction,
    WeightFunctionB::AbstractWeightFunction, BackgroundQuantities::BackgroundQuantities,
    ::AbstractCosmology, CosmologicalGrid::AbstractCosmologicalGrid,
    PowerSpectrum::AbstractPowerSpectrum, ::NumericalIntegrationSimpson)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for idx_a in 1:size(Cℓ.CℓArray, 2)
        for idx_b in idx_a:size(Cℓ.CℓArray, 3)
            for idx_l in 1:size(Cℓ.CℓArray, 1)
                Integrand = c_0 .* WeightFunctionA.WeightFunctionArray[idx_a, :] .*
                WeightFunctionB.WeightFunctionArray[idx_b, :] ./
                (BackgroundQuantities.HZArray .* BackgroundQuantities.χZArray.^2) .*
                PowerSpectrum.InterpolatedPowerSpectrum[idx_l,:]
                Cℓ.CℓArray[idx_l, idx_a, idx_b] = 
                NumericalIntegration.integrate(CosmologicalGrid.ZArray,
                Integrand, SimpsonEvenFast())
                Cℓ.CℓArray[idx_l, idx_b, idx_a] = Cℓ.CℓArray[idx_l, idx_a, idx_b]
            end
        end
    end
end

"""
    ComputeCℓ!(Cℓ::AbstractCℓ, WeightFunctionA::AbstractWeightFunction,
    WeightFunctionB::AbstractWeightFunction,
    BackgroundQuantities::AbstractBackgroundQuantities,
    ::AbstractCosmology, CosmologicalGrid::AbstractCosmologicalGrid,
    PowerSpectrum::AbstractPowerSpectrum, ::CustomSimpson)

This function evaluates the Angular Coefficients for all tomographic bins and
multipole values. In order to evaluate the numerical integrals, it has been
implemented the Simpson rule. The computation is accelerated by
[LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) .
"""
function  ComputeCℓ!(Cℓ::AbstractCℓ, WeightFunctionA::AbstractWeightFunction,
    WeightFunctionB::AbstractWeightFunction,
    BackgroundQuantities::AbstractBackgroundQuantities,
    ::AbstractCosmology, CosmologicalGrid::AbstractCosmologicalGrid,
    PowerSpectrum::AbstractPowerSpectrum, ::CustomSimpson)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    check = true
    while check == true
        SimpsonWeights = SimpsonWeightArray(length(CosmologicalGrid.ZArray))
        ZStep = (last(CosmologicalGrid.ZArray)-first(CosmologicalGrid.ZArray)) /
        (length(CosmologicalGrid.ZArray)-1)
        Cℓ.CℓArray = CustomCℓIntegrator!(SimpsonWeights, Cℓ.CℓArray,
        WeightFunctionA.WeightFunctionArray, WeightFunctionB.WeightFunctionArray,
        ZStep, BackgroundQuantities.HZArray, BackgroundQuantities.χZArray,
        PowerSpectrum.InterpolatedPowerSpectrum)
        if any(isnan,Cℓ.CℓArray)
            
        else
            check = false
        end
    end
end


"""
    CustomCℓIntegrator!(SimpsonWeights::AbstractArray{T}, CℓArray::AbstractArray{T, N},
    WArrayA::AbstractArray{T,N}, WArrayB::AbstractArray{T,N}, ZStep::T,
    HZArray::AbstractArray{T}, χZArray::AbstractArray{T},
    InterpolatedPmm::AbstractArray{T,N}) where {T, N}

This function computes the Cℓ, given the appropriate input, using the Simpson integration 
rule. The computation is accelerated by
[LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl) .
"""
function CustomCℓIntegrator!(SimpsonWeights::AbstractArray{T}, CℓArray::AbstractArray{T, N},
    WArrayA::AbstractArray{T,2}, WArrayB::AbstractArray{T,2}, ZStep::T,
    HZArray::AbstractArray{T}, χZArray::AbstractArray{T},
    InterpolatedPmm::AbstractArray{<:Number}) where {T, N}
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    @avx for i ∈ axes(CℓArray,2),
        j ∈ axes(CℓArray,3),
        l ∈ axes(CℓArray,1)
        for z ∈ axes(χZArray,1)
            CℓArray[l,i,j] += c_0 * WArrayA[i, z] * WArrayB[j, z] /
            (HZArray[z] * χZArray[z]^2) * InterpolatedPmm[l,z] * SimpsonWeights[z]
        end
    end
    CℓArray .*= ZStep
end
    

function  ComputeCℓ!(Cℓ::AbstractCℓ, TransferFunctionA::AbstractTransferFunction,
    TransferFunctionB::AbstractTransferFunction, ::AbstractCosmology,
    CosmologicalGrid::AbstractCosmologicalGrid, PowerSpectrum::AbstractPowerSpectrum, 
    ::BeyondLimber)
    Integrand = zeros(size(Cℓ.CℓArray))
    WeightsMatrix = UnevenTrapzWeightMatrix(CosmologicalGrid.KBeyondLimberArray)
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
