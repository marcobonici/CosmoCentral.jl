function EvaluateFisherMatrixElement!(fishermatrix::Fisherαβ, Cov::aₗₘCovariance,
    ∂Cℓα::∂Cℓ, ∂Cℓβ::∂Cℓ, cosmogrid::CosmologicalGrid, Parα::String, Parβ::String)
    fisherelement = 0
    for ℓidx in 1:length(cosmogrid.MultipolesArray)
        for i in 1:length(∂Cℓα.∂CℓArray[1,:,1])
            for j in 1:length(∂Cℓα.∂CℓArray[1,1,:])
                for m in 1:length(∂Cℓβ.∂CℓArray[1,:,1])
                    for n in 1:length(∂Cℓβ.∂CℓArray[1,1,:])
                        fisherelement += ∂Cℓα.∂CℓArray[ℓidx, i, j] * Cov.Covariance⁻¹[ℓidx, j, m] *
                        ∂Cℓβ.∂CℓArray[ℓidx, m, n] * Cov.Covariance⁻¹[ℓidx, n, i]
                    end
                end
            end
        end
    end
    fishermatrix.FisherDict[Parα*"_"*Parβ] = fisherelement
end