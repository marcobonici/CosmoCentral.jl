function InstantiateEvaluateCovariance(cℓ::AbstractCℓ, ConvDens::AbstractConvolvedDensity,
    cosmogrid::CosmologicalGrid)
    Cov = aₗₘCovariance()
    Cov.Cℓ = cℓ
    Cov.Noise      = similar(cℓ.CℓArray) .* 0.
    Cov.Covariance = similar(cℓ.CℓArray) .* 0.
    EvaluateNoise!(Cov, ConvDens)
    EvaluateCovariance!(Cov, cosmogrid)
    InvertCovariance!(Cov)
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

function EvaluateCovariance!(Cov::aₗₘCovariance, cosmogrid::CosmologicalGrid)
    Cov.Covariance = (Cov.Cℓ.CℓArray + Cov.Noise)
    f_sky = 0.36
    Δℓ =ΔLogSpaced(first(cosmogrid.MultipolesArray), last(cosmogrid.MultipolesArray),
    length(cosmogrid.MultipolesArray))
    for (ℓidx, ℓvalue) in enumerate(cosmogrid.MultipolesArray)
        Cov.Covariance[ℓidx,:,:] *= sqrt(2/((2* ℓvalue+1)*f_sky*Δℓ[ℓidx]))
    end
end

function InvertCovariance!(Cov::aₗₘCovariance)
    for ℓidx in 1:length(Cov.Covariance[:,1,1])
        Cov.Covariance⁻¹[ℓidx,:,:] = inv(Cov.Covariance[ℓidx,:,:])
    end
end