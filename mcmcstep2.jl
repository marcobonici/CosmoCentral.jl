# Load Turing and MCMCChains.
using Turing
using LinearAlgebra
using CosmoCentral

w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()
println("Loaded central cosmology")

@model function gdemo(vecpCℓData, Cov, ConvolvedDensity, EuclidBias, EuclidIA,
    CosmologicalGrid)
    w₀ ~ Uniform(-1.5, -0.5)

    Cosmology = CosmoCentral.Flatw0waCDMCosmology()
    w0waCDMCosmology.w0 = w₀
    CℓMCMC = CosmoCentral.EvaluateCℓMCMCStep(w0waCDMCosmology, ConvolvedDensity, EuclidBias,
    EuclidIA, CosmologicalGrid)

    ℓnumber = length(CℓMCMC.CℓArray[:,1,1])
    inumber = length(CℓMCMC.CℓArray[1,:,1])
    L = CosmoCentral.EliminationMatrix(inumber)

    vecpCℓMCMC = zeros(1,floor(Int,inumber*0.5*(inumber+1)))
    for ℓ in 1:ℓnumber
        vecCℓMCMC = vec(CℓMCMC.CℓArray[ℓ,:,:])
        vecpCℓMCMC = zeros(floor(Int,inumber*0.5*(inumber+1)))
        LinearAlgebra.mul!(vecpCℓMCMC, L, vecCℓMCMC)

        vecpCℓData[ℓ,:] ~ MvNormal(vecpCℓMCMC, Cov.Covariance[ℓ,:,:])
    end
end
MultipolesArrayTemp = CosmoCentral.LogSpaced(10.,3000., 101)
MultipolesArray = zeros(100)
#MultipolesWidths = vcat(CosmoCentral.Difference(MultipolesArrayTemp), ones(2000))
MultipolesWidths = CosmoCentral.Difference(MultipolesArrayTemp)
for i in 1:100
    MultipolesArray[i] = (MultipolesArrayTemp[i+1]+MultipolesArrayTemp[i])/2
end
CosmologicalGrid = CosmoCentral.CosmologicalGrid(ZArray = LinRange(0.001, 4., 500),
KArray = CosmoCentral.LogSpaced(1e-5, 50., 1000), ℓBinCenters = MultipolesArray,
ℓBinWidths = MultipolesWidths);

AnalitycalDensity = CosmoCentral.AnalitycalDensity()
CosmoCentral.NormalizeAnalitycalDensity!(AnalitycalDensity)

InstrumentResponse = CosmoCentral.InstrumentResponse()

ConvolvedDensity = CosmoCentral.ConvolvedDensity(DensityGridArray = ones(10, length(CosmologicalGrid.ZArray)))
CosmoCentral.NormalizeConvolvedDensity!(ConvolvedDensity, AnalitycalDensity, InstrumentResponse, CosmologicalGrid)

CosmoCentral.ComputeConvolvedDensityGrid!(CosmologicalGrid, ConvolvedDensity, AnalitycalDensity, InstrumentResponse)

EuclidBias = CosmoCentral.EuclidBias()
EuclidIA = CosmoCentral.ExtendedNLIA()

CℓData = CosmoCentral.EvaluateCℓMCMCStep(w0waCDMCosmology, ConvolvedDensity, EuclidBias,
EuclidIA, CosmologicalGrid)
#TODO now only LL, but this need to be more flexible...maybe list with probes?
Covaₗₘ = CosmoCentral.InstantiateEvaluateCovariance(CℓData, ConvolvedDensity, CosmologicalGrid, "Lensing",
"Lensing")
CovCℓ = CosmoCentral.InstantiateEvaluateCovariance(Covaₗₘ)
inumber = 10
L = CosmoCentral.EliminationMatrix(inumber)
vecpCℓData = zeros(100,floor(Int,inumber*0.5*(inumber+1)))
for ℓ in 1:length(CℓData.CℓArray[:,1,1])
    vecCℓData = vec(CℓData.CℓArray[ℓ,:,:])
    vecpCℓDataT = zeros(floor(Int,inumber*0.5*(inumber+1)))
    LinearAlgebra.mul!(vecpCℓDataT, L, vecCℓData)
    vecpCℓData[ℓ,:] = vecpCℓDataT
    Matrix = (CovCℓ.Covariance[ℓ,:,:] .+ transpose(CovCℓ.Covariance[ℓ,:,:])) ./2
    CovCℓ.Covariance[ℓ,:,:] = Matrix
end
println(size(CovCℓ.Covariance[1,:,:]))




model = gdemo(vecpCℓData, CovCℓ, ConvolvedDensity, EuclidBias, EuclidIA,
CosmologicalGrid)

chains = sample(model, MH(), 100; save_state = true)
write("first_chain-file.jls", chains)
