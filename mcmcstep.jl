# Load Turing and MCMCChains.
using Turing
using CosmoCentral
using Distributed
using Dates
using ClusterManagers
using PyCall
using LinearAlgebra

n_process = 6
#add  n_processes on the long queue
for i in 1:n_process
    time = string(Dates.now())
    ClusterManagers.addprocs_lsf(1; bsub_flags = `-q long -o Turing/outfiles/$time.out -e Turing/outfiles/$time.err -M 3000 -m g2farm7.ge.infn.it`)
end

@everywhere using Distributed, Turing, LinearAlgebra, PyCall
@everywhere using CosmoCentral

@everywhere w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()
println("Loaded central cosmology")

@everywhere @model function gdemo(vecpCℓData, Cov, ConvolvedDensity, EuclidBias, EuclidIA,
    CosmologicalGrid)
    w₀ ~ Uniform(-1.5, -0.5)
    wₐ ~ Uniform(-0.5, 0.5)
    ΩM ~ Uniform(0.2, 0.4)
    ns ~ Uniform(0.8, 1.1)

    w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()
    w0waCDMCosmology.w0 = w₀
    w0waCDMCosmology.wa = wₐ
    w0waCDMCosmology.ΩM = ΩM
    w0waCDMCosmology.ns = ns
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

@everywhere MultipolesArrayTemp = CosmoCentral.LogSpaced(10.,3000., 101)
@everywhere MultipolesArray = zeros(100)
#MultipolesWidths = vcat(CosmoCentral.Difference(MultipolesArrayTemp), ones(2000))
@everywhere MultipolesWidths = CosmoCentral.Difference(MultipolesArrayTemp)
@everywhere for i in 1:100
    MultipolesArray[i] = (MultipolesArrayTemp[i+1]+MultipolesArrayTemp[i])/2
end
@everywhere CosmologicalGrid = CosmoCentral.CosmologicalGrid(ZArray = LinRange(0.001, 4., 500),
KArray = CosmoCentral.LogSpaced(1e-5, 50., 1000), ℓBinCenters = MultipolesArray,
ℓBinWidths = MultipolesWidths);

@everywhere AnalitycalDensity = CosmoCentral.AnalitycalDensity()
@everywhere CosmoCentral.NormalizeAnalitycalDensity!(AnalitycalDensity)

@everywhere InstrumentResponse = CosmoCentral.InstrumentResponse()

@everywhere ConvolvedDensity = CosmoCentral.ConvolvedDensity(DensityGridArray = ones(10, length(CosmologicalGrid.ZArray)))
@everywhere CosmoCentral.NormalizeConvolvedDensity!(ConvolvedDensity, AnalitycalDensity, InstrumentResponse, CosmologicalGrid)

@everywhere CosmoCentral.ComputeConvolvedDensityGrid!(CosmologicalGrid, ConvolvedDensity, AnalitycalDensity, InstrumentResponse)

@everywhere EuclidBias = CosmoCentral.EuclidBias()
@everywhere EuclidIA = CosmoCentral.ExtendedNLIA()

@everywhere CℓData = CosmoCentral.EvaluateCℓMCMCStep(w0waCDMCosmology, ConvolvedDensity, EuclidBias,
EuclidIA, CosmologicalGrid)
#TODO now only LL, but this need to be more flexible...maybe list with probes?
@everywhere Covaₗₘ = CosmoCentral.InstantiateEvaluateCovariance(CℓData, ConvolvedDensity, CosmologicalGrid, "Lensing",
"Lensing")
@everywhere CovCℓ = CosmoCentral.InstantiateEvaluateCovariance(Covaₗₘ)
@everywhere inumber = 10
@everywhere L = CosmoCentral.EliminationMatrix(inumber)
@everywhere vecpCℓData = zeros(100,floor(Int,inumber*0.5*(inumber+1)))
@everywhere for ℓ in 1:length(CℓData.CℓArray[:,1,1])
    vecCℓData = vec(CℓData.CℓArray[ℓ,:,:])
    vecpCℓDataT = zeros(floor(Int,inumber*0.5*(inumber+1)))
    LinearAlgebra.mul!(vecpCℓDataT, L, vecCℓData)
    vecpCℓData[ℓ,:] = vecpCℓDataT
    Matrix = (CovCℓ.Covariance[ℓ,:,:] .+ transpose(CovCℓ.Covariance[ℓ,:,:])) ./2
    CovCℓ.Covariance[ℓ,:,:] = Matrix
end
println(size(CovCℓ.Covariance[1,:,:]))




@everywhere model = gdemo(vecpCℓData, CovCℓ, ConvolvedDensity, EuclidBias, EuclidIA,
CosmologicalGrid)

chains = sample(model, MH(), MCMCDistributed(), 10, n_process; save_state = true)
write("first_chain-file.jls", chains)