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
function  ComputeCℓ!(Cℓ::AbstractCℓ, WFA::AbstractWeightFunction,
    WFB::AbstractWeightFunction, backquant::BackgroundQuantities,
    ::AbstractCosmology, cosmogrid::AbstractCosmologicalGrid,
    pmm::AbstractPowerSpectrum, ::NumericalIntegrationSimpson)
    c_0 = 2.99792458e5 #TODO: find a package containing the exact value of
                       #physical constants involved in calculations
    for idx_a in 1:size(Cℓ.CℓArray, 2)
        for idx_b in idx_a:size(Cℓ.CℓArray, 3)
            for idx_l in 1:size(Cℓ.CℓArray, 1)
                Integrand = c_0 .* WFA.WeightFunctionArray[idx_a, :] .*
                WFB.WeightFunctionArray[idx_b, :] ./
                (backquant.HZArray .* backquant.χZArray.^2) .*
                pmm.InterpolatedPowerSpectrum[idx_l,:]
                Cℓ.CℓArray[idx_l, idx_a, idx_b] = 
                NumericalIntegration.integrate(cosmogrid.ZArray,
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
function  ComputeCℓ!(Cℓ::AbstractCℓ, WFA::AbstractWeightFunction,
    WFB::AbstractWeightFunction,
    backquant::AbstractBackgroundQuantities,
    ::AbstractCosmology, cosmogrid::AbstractCosmologicalGrid,
    pmm::AbstractPowerSpectrum, ::CustomSimpson)
    c_0 = 2.99792458e5 #TODO: #108 find a package containing the exact value of
                       #physical constants involved in calculations
    check = true
    while check == true
        SimpsonWeights = SimpsonWeightArray(length(cosmogrid.ZArray))
        ZStep = (last(cosmogrid.ZArray)-first(cosmogrid.ZArray)) /
        (length(cosmogrid.ZArray)-1)
        Cℓ.CℓArray = CustomCℓIntegrator!(SimpsonWeights, Cℓ.CℓArray,
        WFA.WeightFunctionArray, WFB.WeightFunctionArray,
        ZStep, backquant.HZArray, backquant.χZArray,
        pmm.InterpolatedPowerSpectrum)
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

This function computes the Cℓ for two probes A and B using the Simpson integration 
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
    

function  ComputeCℓ!(Cℓ::AbstractCℓ, TFA::AbstractTransferFunction,
    TFB::AbstractTransferFunction, ::AbstractCosmology,
    cosmogrid::AbstractCosmologicalGrid, pmm::AbstractPowerSpectrum, 
    ::BeyondLimber)
    Integrand = zeros(size(Cℓ.CℓArray))
    WeightsMatrix = UnevenTrapzWeightMatrix(cosmogrid.KBeyondLimberArray)
    @avx for i ∈ axes(Cℓ.CℓArray,2),
        j ∈ axes(Cℓ.CℓArray,3),
        l ∈ axes(Cℓ.CℓArray,1)
        for k ∈ axes(cosmogrid.KBeyondLimberArray,2)
            Integrand[l,i,j] += 2 / π * WeightsMatrix[l, k] *
            TFA.TransferFunctionArray[i, l, k] *
            TFB.TransferFunctionArray[j, l, k] *
            pmm.InterpolatedPowerSpectrumBeyondLimber[l,k] *
            cosmogrid.KBeyondLimberArray[l, k]^2
        end
    end
    Cℓ.CℓArray = Integrand
end
