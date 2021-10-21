function EvaluateFisherMatrixElement!(FisherMatrix::Fisherαβ, Cov::aₗₘCovariance,
    ∂Cℓα::∂Cℓ, ∂Cℓβ::∂Cℓ, Parα::String, Parβ::String)
    αMatrix = mymatmul(∂Cℓα.∂CℓArray,Cov.Covariance⁻¹)
    βMatrix = mymatmul(∂Cℓβ.∂CℓArray,Cov.Covariance⁻¹)
    Fisherℓ = mymatmul(αMatrix, βMatrix)
    fisherelementℓ, fisherelement = SumℓAndTrace(Fisherℓ)
    FisherMatrix.FisherDict[Parα*"_"*Parβ] = fisherelement
    FisherMatrix.FisherℓDict[Parα*"_"*Parβ] = fisherelementℓ
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

function SumℓAndTrace(Fisherℓ::Array{Float64, 3})
    fisherelement = 0
    fisherelementℓ = zeros(size(Fisherℓ)[1])
    @avx for ℓ in 1:size(Fisherℓ)[1]
        for i in 1:size(Fisherℓ)[2]
            fisherelementℓ[ℓ] += Fisherℓ[ℓ, i, i]
        end
    end
    return fisherelementℓ, sum(fisherelementℓ)
end

function EvaluateFisherMatrixElement!(FisherMatrix::Fisherαβ, Cov::CℓCovariance,
    ∂Cℓα::∂Cℓ, ∂Cℓβ::∂Cℓ, Parα::String, Parβ::String)
    ℓnumber = length(∂Cℓα.∂CℓArray[:,1,1])
    inumber = length(∂Cℓα.∂CℓArray[1,:,1])
    L = EliminationMatrix(inumber)

    vecp∂CℓαᵀCov⁻¹ = zeros(1,floor(Int,inumber*0.5*(inumber+1)))
    fisher_temp = zeros(1,1)
    FisherMatrix.FisherℓDict[Parα*"_"*Parβ] = zeros(ℓnumber)
    fisherelement = 0
    for ℓ in 1:ℓnumber
        vec∂Cℓα = vec(∂Cℓα.∂CℓArray[ℓ,:,:])
        vecp∂Cℓα = zeros(floor(Int,inumber*0.5*(inumber+1)))
        LinearAlgebra.mul!(vecp∂Cℓα, L, vec∂Cℓα)
        
        vec∂Cℓβ = vec(∂Cℓβ.∂CℓArray[ℓ,:,:])
        vecp∂Cℓβ = zeros(floor(Int,inumber*0.5*(inumber+1)))
        LinearAlgebra.mul!(vecp∂Cℓβ, L, vec∂Cℓβ)

        Covariance⁻¹ = Cov.Covariance⁻¹[ℓ,:,:]
        LinearAlgebra.mul!(vecp∂CℓαᵀCov⁻¹, transpose(vecp∂Cℓα), Covariance⁻¹)
        LinearAlgebra.mul!(fisher_temp, vecp∂CℓαᵀCov⁻¹, vecp∂Cℓβ)
        FisherMatrix.FisherℓDict[Parα*"_"*Parβ][ℓ] = fisher_temp[1,1]
        fisherelement += fisher_temp[1,1]
    end
    FisherMatrix.FisherDict[Parα*"_"*Parβ] = fisherelement
end