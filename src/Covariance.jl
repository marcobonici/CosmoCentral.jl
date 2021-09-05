"""
    InstantiateEvaluateCovariance(cℓ::AbstractCℓ, ConvDens::AbstractConvolvedDensity,
    cosmogrid::CosmologicalGrid, ProbeA::String, ProbeB::String)

This function evaluates and returns the [`aₗₘCovariance`](@ref), according to the following
formula:
```math
\\Sigma_{i j}^{\\mathrm{AB}}(\\ell)=\\sqrt{\\frac{2}{(2 \\ell+1) \\Delta \\ell
f_{\\mathrm{sky}}}}\\left(C_{i j}^{\\mathrm{AB}}(\\ell)+N_{i j}^{\\mathrm{AB}}(\\ell)
\\right),
```
where ``\\mathrm{A}`` and ``\\mathrm{B}`` are the probes, ``\\mathrm{i}`` and
``\\mathrm{j}`` are the tomographic bins, ``\\ell`` and ``\\Delta\\ell`` are respectively the
central value and width of the considered multipole bin, ``f_{\\mathrm{sky}}`` is the sky
fraction covered by the field considered, ``N_{i j}^{\\mathrm{AB}}`` is the noise matrix.
"""
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

"""
    InstantiateEvaluateCovariance(Covaₗₘ::aₗₘCovariance)

This function evaluates and returns the [`CℓCovariance`](@ref), using the following formula:



```math
\\boldsymbol{\\Xi}(\\ell)= \\left(\\boldsymbol{D}_{n}^{T}\\left( \\left(\\boldsymbol{\\Sigma}(\\ell)\\right)^{-1} \\otimes
\\left(\\boldsymbol{\\Sigma}(\\ell)\\right)^{-1}\\right) \\boldsymbol{D}_{n}\\right)^{-1}
```
where ``\\boldsymbol{D}_{n}`` is the [`DuplicationMatrix`](@ref).
"""
function InstantiateEvaluateCovariance(Covaₗₘ::aₗₘCovariance)
    Cov = CℓCovariance()
    ℓnumber = length(Covaₗₘ.Cℓ.CℓArray[:,1,1])
    inumber = length(Covaₗₘ.Cℓ.CℓArray[1,1,:])
    D = DuplicationMatrix(inumber)
    Dᵀ = Transpose(D)
    kronCovaₗₘ = zeros(inumber^2,inumber^2)
    DᵀkronCovaₗₘ = zeros(floor(Int,inumber*0.5*(inumber+1)),inumber^2)
    Cov.Covariance = zeros(ℓnumber, floor(Int,inumber*0.5*(inumber+1)),
    floor(Int,inumber*0.5*(inumber+1)))
    Cov.Covariance = zeros(size(Cov.Covariance))
    Cov.Covariance⁻¹ = zeros(size(Cov.Covariance))
    TempCov = zeros(floor(Int,inumber*0.5*(inumber+1)),
    floor(Int,inumber*0.5*(inumber+1)))
    for ℓ in 1:ℓnumber
        kronCovaₗₘ = kron(Covaₗₘ.Covariance⁻¹[ℓ,:,:], Covaₗₘ.Covariance⁻¹[ℓ,:,:])
        LinearAlgebra.mul!(DᵀkronCovaₗₘ, Dᵀ, kronCovaₗₘ)
        LinearAlgebra.mul!(TempCov, DᵀkronCovaₗₘ, D)
        Cov.Covariance⁻¹[ℓ,:,:] = TempCov 
        Cov.Covariance[ℓ,:,:] = inv(Cov.Covariance⁻¹[ℓ,:,:])
        Matrix = (Cov.Covariance[ℓ,:,:] .+ transpose(Cov.Covariance[ℓ,:,:])) ./2
        Cov.Covariance[ℓ,:,:] = Matrix
        Cov.Covariance⁻¹[ℓ,:,:] = inv(Cov.Covariance[ℓ,:,:])
    end
    return Cov
end