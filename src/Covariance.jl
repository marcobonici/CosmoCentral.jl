function InstantiateEvaluateCovariance(cℓ::AbstractCℓ, ConvDens::AbstractConvolvedDensity)
    Cov = aₗₘCovariance()
    Cov.Cℓ = cℓ
    Cov.Noise      = similar(cℓ.CℓArray) .* 0.
    Cov.Covariance = similar(cℓ.CℓArray) .* 0.
    EvaluateNoise!(Cov, ConvDens)
    EvaluateCovariance!(Cov)
    return Cov
end

function EvaluateNoise!(Cov::aₗₘCovariance, ConvDens::AbstractConvolvedDensity)
    for lidx in 1:length(Cov.Cℓ.CℓArray[:,1,1])
        for iidx in 1:length(Cov.Cℓ.CℓArray[1,:,1])
            Cov.Noise[lidx, iidx, iidx] = 1 ./ 
            (ConvDens.SurfaceDensityArray[iidx]*3437.746771^2)
        end
    end
    #TODO: the multiplication here is to convert from square degrees to steradians. This
    #should be done better...
end

function EvaluateCovariance!(Cov::aₗₘCovariance)
    Cov.Covariance = Cov.Cℓ.CℓArray + Cov.Noise
end