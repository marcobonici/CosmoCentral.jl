# Load Turing and MCMCChains.
using Turing
using LinearAlgebra
using CosmoCentral
using BSON
using Flux

w0waCDMCosmology = CosmoCentral.w0waCDMCosmology()
println("Loaded central cosmology")

cosmo = BSON.load("cosmoemulatorll.bson", @__MODULE__)
CosmoEmulatorLL = cosmo[:cosmo]
cosmo = BSON.load("cosmoemulatorgg.bson", @__MODULE__)
CosmoEmulatorGG = cosmo[:cosmo]
println("Loaded Emulators")

ProbeA = "Lensing"

ProbeB = "PhotometricClustering"

@model function gdemo(vecpCℓDataLL, CovLL, CosmoEmulatorLL, vecpCℓDataGG, CovGG,
    CosmoEmulatorGG, CosmologicalGrid)
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

    CℓMCMCLL = CosmoCentral.Cℓ(CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters), 10,10))
    CℓMCMCLL.CℓArray = CosmoCentral.ComputeCℓ(CosmoEmulatorLL, w0waCDMCosmology, CosmologicalGrid)

    CℓMCMCGG = CosmoCentral.Cℓ(CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters), 10,10))
    CℓMCMCGG.CℓArray = CosmoCentral.ComputeCℓ(CosmoEmulatorGG, w0waCDMCosmology, CosmologicalGrid)

    ℓnumber = length(CℓMCMCLL.CℓArray[:,1,1])
    inumber = length(CℓMCMCLL.CℓArray[1,:,1])
    #L = CosmoCentral.EliminationMatrix(inumber)
    vecpCℓMCMCLL = zeros(1,floor(Int,inumber*0.5*(inumber+1)))
    vecpCℓMCMCGG = zeros(1,floor(Int,inumber*0.5*(inumber+1)))

    for ℓ in 1:ℓnumber
        vecCℓMCMCLL = vec(CℓMCMCLL.CℓArray[ℓ,:,:])
        vecpCℓMCMCLL = zeros(floor(Int,inumber*0.5*(inumber+1)))
        LinearAlgebra.mul!(vecpCℓMCMCLL, L, vecCℓMCMCLL)

        vecpCℓDataLL[ℓ,:] ~ MvNormal(vecpCℓMCMCLL, CovLL.Covariance[ℓ,:,:])

        vecCℓMCMCGG = vec(CℓMCMCGG.CℓArray[ℓ,:,:])
        vecpCℓMCMCGG = zeros(floor(Int,inumber*0.5*(inumber+1)))
        LinearAlgebra.mul!(vecpCℓMCMCGG, L, vecCℓMCMCGG)

        vecpCℓDataGG[ℓ,:] ~ MvNormal(vecpCℓMCMCGG, CovGG.Covariance[ℓ,:,:])
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
#CosmoCentral.NormalizeConvolvedDensity!(ConvolvedDensity, AnalitycalDensity, InstrumentResponse, CosmologicalGrid)

#CosmoCentral.ComputeConvolvedDensityGrid!(CosmologicalGrid, ConvolvedDensity, AnalitycalDensity, InstrumentResponse)

#EuclidBias = CosmoCentral.EuclidBias()
#EuclidIA = CosmoCentral.ExtendedNLIA()

CℓDataLL = CosmoCentral.Cℓ(CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters), 10,10))
CℓDataLL.CℓArray = CosmoCentral.ComputeCℓ(CosmoEmulatorLL, w0waCDMCosmology, CosmologicalGrid)
CovaₗₘLL = CosmoCentral.InstantiateEvaluateCovariance(CℓDataLL, ConvolvedDensity,
CosmologicalGrid, ProbeA, ProbeA)

CℓDataGG = CosmoCentral.Cℓ(CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters), 10,10))
CℓDataGG.CℓArray = CosmoCentral.ComputeCℓ(CosmoEmulatorGG, w0waCDMCosmology, CosmologicalGrid)
CovaₗₘGG = CosmoCentral.InstantiateEvaluateCovariance(CℓDataGG, ConvolvedDensity,
CosmologicalGrid, ProbeB, ProbeB)



println("Created basic covariances")
CovCℓLL = CosmoCentral.InstantiateEvaluateCovariance(CovaₗₘLL)
inumber = 10
const L = CosmoCentral.EliminationMatrix(inumber)
vecpCℓDataLL = zeros(30,floor(Int,inumber*0.5*(inumber+1)))
for ℓ in 1:length(CℓDataLL.CℓArray[:,1,1])
    vecCℓDataLL = vec(CℓDataLL.CℓArray[ℓ,:,:])
    vecpCℓDataLLT = zeros(floor(Int,inumber*0.5*(inumber+1)))
    LinearAlgebra.mul!(vecpCℓDataLLT, L, vecCℓDataLL)
    vecpCℓDataLL[ℓ,:] = vecpCℓDataLLT
    Matrix = (CovCℓLL.Covariance[ℓ,:,:] .+ transpose(CovCℓLL.Covariance[ℓ,:,:])) ./2
    CovCℓLL.Covariance[ℓ,:,:] = Matrix
end

CovCℓGG = CosmoCentral.InstantiateEvaluateCovariance(CovaₗₘGG)
inumber = 10
const L = CosmoCentral.EliminationMatrix(inumber)
vecpCℓDataGG = zeros(30,floor(Int,inumber*0.5*(inumber+1)))
for ℓ in 1:length(CℓDataGG.CℓArray[:,1,1])
    vecCℓDataGG = vec(CℓDataGG.CℓArray[ℓ,:,:])
    vecpCℓDataGGT = zeros(floor(Int,inumber*0.5*(inumber+1)))
    LinearAlgebra.mul!(vecpCℓDataGGT, L, vecCℓDataGG)
    vecpCℓDataGG[ℓ,:] = vecpCℓDataGGT
    Matrix = (CovCℓGG.Covariance[ℓ,:,:] .+ transpose(CovCℓGG.Covariance[ℓ,:,:])) ./2
    CovCℓGG.Covariance[ℓ,:,:] = Matrix
end

println(size(CovCℓLL.Covariance))
println(size(vecpCℓDataLL))
println("Created Covariance Matrix and mock data")



model = gdemo(vecpCℓDataLL, CovCℓLL, CosmoEmulatorLL, vecpCℓDataGG, CovCℓGG, CosmoEmulatorGG,
CosmologicalGrid)

chains = sample(model, MH(), MCMCThreads(), 1000, 100; save_state = false)
write("first_chain-file.jls", chains)
