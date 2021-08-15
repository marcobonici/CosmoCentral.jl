function InstantiateEvaluateCovariance(cℓ::AbstractCℓ, ConvDens::AbstractConvolvedDensity,
    cosmogrid::CosmologicalGrid, ProbeA::String, ProbeB::String)
    Cov = aₗₘCovariance()
    Cov.Cℓ = cℓ
    Cov.Noise      = zeros(size(cℓ.CℓArray))
    Cov.Covariance = zeros(size(cℓ.CℓArray))
    EvaluateNoise!(Cov, ConvDens, ProbeA::String, ProbeB::String)
    EvaluateCovariance!(Cov, cosmogrid)
    InvertCovariance!(Cov)
    return Cov
end

function EvaluateNoise!(Cov::aₗₘCovariance, ConvDens::AbstractConvolvedDensity,
    ProbeA::String, ProbeB::String)
    if ProbeA == ProbeB
        if ProbeA =="Lensing"
            probe_factor = 0.3
        else
            probe_factor = 1.0
        end
        for lidx in 1:length(Cov.Cℓ.CℓArray[:,1,1])
            for iidx in 1:length(Cov.Cℓ.CℓArray[1,:,1])
                Cov.Noise[lidx, iidx, iidx] = probe_factor^2 ./ 
                (ConvDens.SurfaceDensityArray[iidx]*3437.746771^2)
                #TODO: the multiplication here is to convert from square degrees to
                #steradians. This should be done better...
            end
        end
    end
end

function EvaluateCovariance!(Cov::aₗₘCovariance, cosmogrid::CosmologicalGrid)
    Cov.Covariance = (Cov.Cℓ.CℓArray + Cov.Noise)
    f_sky = 0.363610260832152
    #TODO: this is hardcoded...
    for (ℓidx, ℓvalue) in enumerate(cosmogrid.ℓBinCenters)
        Cov.Covariance[ℓidx,:,:] *= sqrt(2/((2* ℓvalue+1)*f_sky*cosmogrid.ℓBinWidths[ℓidx]))
    end
end

function InvertCovariance!(Cov::AbstractCovariance)
    for ℓidx in 1:length(Cov.Covariance[:,1,1])
        Cov.Covariance⁻¹[ℓidx,:,:] = inv(Cov.Covariance[ℓidx,:,:])
    end
end

function InstantiateEvaluateCovariance(Covaₗₘ::aₗₘCovariance)
    Cov = CℓCovariance()
    ℓnumber = length(Covaₗₘ.Cℓ.CℓArray[:,1,1])
    inumber = length(Covaₗₘ.Cℓ.CℓArray[1,1,:])
    D = DuplicationMatrix(inumber)
    Dᵀ = Transpose(D)
    kronCovaₗₘ = zeros(inumber^2,inumber^2)
    DᵀkronCovaₗₘ = zeros( floor(Int,inumber*0.5*(inumber+1)),inumber^2)
    Cov.Covariance = zeros(ℓnumber, floor(Int,inumber*0.5*(inumber+1)),
    floor(Int,inumber*0.5*(inumber+1)))
    Cov.Covariance = zeros(size(Cov.Covariance))
    Cov.Covariance⁻¹ = zeros(size(Cov.Covariance))
    TempCov = zeros(floor(Int,inumber*0.5*(inumber+1)),
    floor(Int,inumber*0.5*(inumber+1)))
    for ℓ in 1:ℓnumber
        kronCovaₗₘ = kron(Covaₗₘ.Covariance[ℓ,:,:], Covaₗₘ.Covariance[ℓ,:,:])
        LinearAlgebra.mul!(DᵀkronCovaₗₘ, Dᵀ, kronCovaₗₘ)
        LinearAlgebra.mul!(TempCov, DᵀkronCovaₗₘ, D)
        Cov.Covariance[ℓ,:,:] = TempCov
        Cov.Covariance⁻¹[ℓ,:,:] = inv(Cov.Covariance[ℓ,:,:])
    end
    return Cov
end