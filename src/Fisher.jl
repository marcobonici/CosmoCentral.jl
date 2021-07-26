function EvaluateFisherMatrixElement!(fishermatrix::Fisherαβ, Cov::aₗₘCovariance,
    ∂Cℓα::∂Cℓ, ∂Cℓβ::∂Cℓ, cosmogrid::CosmologicalGrid, Parα::String, Parβ::String)
    AMatrix = mymatmul(∂Cℓα.∂CℓArray,Cov.Covariance⁻¹)
    BMatrix = mymatmul(∂Cℓβ.∂CℓArray,Cov.Covariance⁻¹)
    Fisherℓ = mymatmul(AMatrix, BMatrix)
    fisherelement = 0
    @avx for ℓ in 1:size(Fisherℓ)[1]
        for i in 1:size(Fisherℓ)[2]
            fisherelement += Fisherℓ[ℓ, i, i]
        end
    end
    return fisherelement
end

function mymatmul(A::Array{Float64, 3}, B::Array{Float64, 3})
    C = zeros(size(A))
    @avx for l in 1:length(A[:,1,1])
        for i in 1:length(A[1,:,1])
            for j in 1:length(B[1,1,:])
                for k in 1:length(B[1,:,1])
                    C[l,i,j] += A[l,i,k]*B[l,k,j]
                end
            end
        end
    end
    return C
end