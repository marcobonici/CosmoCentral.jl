function EvaluateFisherMatrixElement!(FisherMatrix::Fisherαβ, Cov::aₗₘCovariance,
    ∂Cℓα::∂Cℓ, ∂Cℓβ::∂Cℓ, Parα::String, Parβ::String)
    αMatrix = mymatmul(∂Cℓα.∂CℓArray,Cov.Covariance⁻¹)
    βMatrix = mymatmul(∂Cℓβ.∂CℓArray,Cov.Covariance⁻¹)
    Fisherℓ = mymatmul(αMatrix, βMatrix)
    fisherelementℓ, fisherelement = SumℓAndTrace(Fisherℓ)
    FisherMatrix.FisherDict[Parα*"_"*Parβ] = fisherelement
    FisherMatrix.FisherℓDict[Parα*"_"*Parβ] = fisherelementℓ
end

function mymatmul(A::AbstractArray{T, 3}, B::AbstractArray{T, 3}) where T
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

function SumℓAndTrace(Fisherℓ::AbstractArray{T, 3}) where T
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

#######################
#The code here need a bit of polishing and generalizations. However, it works fine for the
#cases actually tested
#######################

function SumFisher(FisherA::Fisherαβ, FisherB::Fisherαβ)
    FisherC = Fisherαβ()
    if FisherA.ParametersList != FisherB.ParametersList
        error("The two parameters lists are different.")
        #TODO: this can be generalized
    end
    FisherC.FisherMatrix = zeros(length(FisherA.ParametersList),
    length(FisherA.ParametersList))
    FisherC.ParametersList = FisherA.ParametersList
    FisherC.SelectedParametersList = union(FisherA.SelectedParametersList,
    FisherB.SelectedParametersList)
    for (key, value) in FisherA.FisherℓDict
        FisherC.FisherDict[key] = FisherA.FisherDict[key] + FisherB.FisherDict[key]
    end
    for (indexα, Parα) in enumerate(FisherC.SelectedParametersList)
        for (indexβ, Parβ) in enumerate(FisherC.SelectedParametersList)
            FisherC.FisherMatrix[indexα, indexβ] = FisherC.FisherDict[Parα*"_"*Parβ]
        end
    end
    """
    if size(FisherA.FisherMatrixCumℓ[1]) == size(FisherB.FisherMatrixCumℓ[1] != 1)
        FisherC.FisherMatrixCumℓ = zeros(size(FisherA.FisherMatrixCumℓ)[1],
            length(FisherA.ParametersList), length(FisherA.ParametersList) )
        for (key, value) in FisherA.FisherℓDict
            FisherC.FisherℓDict[key] = FisherA.FisherℓDict[key] + FisherB.FisherℓDict[key]
        end
        for (indexα, Parα) in enumerate(FisherC.SelectedParametersList)
            for (indexβ, Parβ) in enumerate(FisherC.SelectedParametersList)
                FisherC.FisherMatrixCumℓ[:, indexα, indexβ] = 
                cumsum(FisherC.FisherℓDict[Parα*"_"*Parβ])
            end
        end
    end
    """
    CosmoCentral.SelectMatrixAndMarginalize!(FisherC.ParametersList, FisherC)
    println(FisherC.ParametersList)
    return FisherC
end

function SelectCorrelationMatrix(Fisher::Fisherαβ, Parameters)
    corr_matrix = zeros(length(Parameters), length(Parameters))
    for (idxα, parα) in enumerate(Parameters)
        for (idxβ, parβ) in enumerate(Parameters)
            corr_matrix[idxα, idxβ] = Fisher.CorrelationMatrixDict[parα*"_"*parβ]
        end
    end
    return corr_matrix
end

function EvaluateFoM(Fisher::Fisherαβ, Parα::String, Parβ::String)
    return sqrt(1/det(SelectCorrelationMatrix(Fisher, [Parα, Parβ])))
end

function RearrangeFisherℓ(Fisher::Fisherαβ, CosmoGrid::CosmologicalGrid, ℓmin, ℓmax)
    FisherNew = deepcopy(Fisher)
    FisherNew.FisherMatrix = zeros(length(Fisher.ParametersList), length(Fisher.ParametersList))
    FisherNew.FisherMatrixCumℓ = zeros(1,1,1)
    for (idxα, parα) in enumerate(Fisher.ParametersList)
        for (idxβ, parβ) in enumerate(Fisher.ParametersList)
            interpFisher = Dierckx.Spline1D(CosmoGrid.ℓBinCenters, Fisher.FisherℓDict[parα*"_"*parβ] ./ CosmoGrid.ℓBinWidths)
            int, err = quadgk(ℓ -> interpFisher(ℓ), ℓmin, ℓmax, rtol=1e-12)
            FisherNew.FisherDict[parα*"_"*parβ] = int
            FisherNew.FisherMatrix[idxα, idxβ] = int
        end
    end
    return FisherNew
end