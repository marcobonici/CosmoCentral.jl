# Load Turing and MCMCChains.
using Turing
using LinearAlgebra
using CosmoCentral
using BSON
using Flux

w0waCDMCosmology = CosmoCentral.w0waCDMCosmology()
println("Loaded central cosmology")

cosmo = BSON.load("cosmoemulatorll.bson", @__MODULE__)
CosmoEmulator = cosmo[:cosmo]
cosmo = BSON.load("cosmoemulatorgg.bson", @__MODULE__)
CosmoEmulatorGG = cosmo[:cosmo]
println("Loaded Emulators")

ProbeA = "PhotometricClustering"

ProbeB = "PhotometricClustering"

@model function gdemo(vecpCℓData, Cov, CosmoEmulator,
    CosmologicalGrid)
    w₀ ~ Uniform(-1.32, -0.72)
    wₐ ~ Uniform(-0.68, 0.72)
    ns ~ Uniform(0.92, 1.0)
    ΩM ~ Uniform(0.29, 0.35)
    ΩB ~ Uniform(0.04, 0.06)  
    σ8 ~ Uniform(0.6,  1.0)
    H0 ~ Uniform(61.,  73.)

    w0waCDMCosmology = CosmoCentral.w0waCDMCosmology()
    w0waCDMCosmology.w0 = w₀
    w0waCDMCosmology.wa = wₐ
    w0waCDMCosmology.ΩM = ΩM
    w0waCDMCosmology.ns = ns
    w0waCDMCosmology.ΩB = ΩB
    w0waCDMCosmology.σ8 = σ8
    w0waCDMCosmology.H0 = H0

    CℓMCMC = CosmoCentral.Cℓ(CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters), 10,10))
    CℓMCMC.CℓArray = CosmoCentral.ComputeCℓ(CosmoEmulator, w0waCDMCosmology, CosmologicalGrid)

    ℓnumber = length(CℓMCMC.CℓArray[:,1,1])
    inumber = length(CℓMCMC.CℓArray[1,:,1])
    #L = CosmoCentral.EliminationMatrix(inumber)
    vecpCℓMCMC = zeros(1,floor(Int,inumber*0.5*(inumber+1)))
    for ℓ in 1:ℓnumber
        vecCℓMCMC = vec(CℓMCMC.CℓArray[ℓ,:,:])
        vecpCℓMCMC = zeros(floor(Int,inumber*0.5*(inumber+1)))
        LinearAlgebra.mul!(vecpCℓMCMC, L, vecCℓMCMC)

        vecpCℓData[ℓ,:] ~ MvNormal(vecpCℓMCMC, Cov.Covariance[ℓ,:,:])
    end
end

MultipolesArrayTemp = CosmoCentral.LogSpaced(10.,3000., 31)
MultipolesArray = zeros(30)
#MultipolesWidths = vcat(CosmoCentral.Difference(MultipolesArrayTemp), ones(2000))
MultipolesWidths = CosmoCentral.Difference(MultipolesArrayTemp)
for i in 1:30
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

#EuclidBias = CosmoCentral.EuclidBias()
#EuclidIA = CosmoCentral.ExtendedNLIA()

CℓData = CosmoCentral.Cℓ(CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters), 10,10))
CℓData.CℓArray = CosmoCentral.ComputeCℓ(CosmoEmulator, w0waCDMCosmology, CosmologicalGrid)
Covaₗₘ = CosmoCentral.InstantiateEvaluateCovariance(CℓData, ConvolvedDensity,
CosmologicalGrid, ProbeA, ProbeB)



println("Created basic covariances")
CovCℓ = CosmoCentral.InstantiateEvaluateCovariance(Covaₗₘ)
inumber = 10
const L = CosmoCentral.EliminationMatrix(inumber)
vecpCℓData = zeros(30,floor(Int,inumber*0.5*(inumber+1)))
for ℓ in 1:length(CℓData.CℓArray[:,1,1])
    vecCℓData = vec(CℓData.CℓArray[ℓ,:,:])
    vecpCℓDataT = zeros(floor(Int,inumber*0.5*(inumber+1)))
    LinearAlgebra.mul!(vecpCℓDataT, L, vecCℓData)
    vecpCℓData[ℓ,:] = vecpCℓDataT
    Matrix = (CovCℓ.Covariance[ℓ,:,:] .+ transpose(CovCℓ.Covariance[ℓ,:,:])) ./2
    CovCℓ.Covariance[ℓ,:,:] = Matrix
end
println(size(CovCℓ.Covariance))
println(size(vecpCℓData))
println("Created Covariance Matrix and mock data")



model = gdemo(vecpCℓData, CovCℓ, CosmoEmulator, CosmologicalGrid)

chains = sample(model, MH(), MCMCThreads(), 1000, 100; save_state = false)
write("first_chain-file.jls", chains)
