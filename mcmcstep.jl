# Load Turing and MCMCChains.
using Turing

using CosmoCentral
using Distributed
using Dates
using ClusterManagers

n_process = 10
#add  n_processes on the long queue
for i in 1:n_process
    time = string(Dates.now())
    ClusterManagers.addprocs_lsf(1; bsub_flags = `-q long -o Turing/outfiles/$time.out -e Turing/outfiles/$time.err -n 2 `)
end



@everywhere using Turing
@everywhere using CosmoCentral
@everywhere w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()
println("Loaded central cosmology")

# Define a model on all processes.
@everywhere @model function gdemo(CℓData, Cov, ConvolvedDensity, EuclidBias, EuclidIA,
    CosmologicalGrid)
    w₀ ~ Uniform(-0.5, -1.5)
    Cosmology = Flatw0waCDMCosmology()
    w0waCDMCosmology.w0 = w₀
    CℓMCMC = CosmoCentral.EvaluateCℓMCMCStep(w0waCDMCosmology, ConvolvedDensity, EuclidBias,
    EuclidIA, CosmologicalGrid)

    ℓnumber = length(CℓData.CℓArray[:,1,1])
    inumber = length(CℓData.CℓArray[1,:,1])
    L = EliminationMatrix(inumber)

    vecpCℓData = zeros(1,floor(Int,inumber*0.5*(inumber+1)))
    vecpCℓMCMC = zeros(1,floor(Int,inumber*0.5*(inumber+1)))
    for ℓ in 1:ℓnumber
        vecCℓData = vec(CℓData.CℓArray[ℓ,:,:])
        vecpCℓData = zeros(floor(Int,inumber*0.5*(inumber+1)))
        LinearAlgebra.mul!(vecpCℓData, L, vecCℓData)

        vecCℓMCMC = vec(CℓMCMC.CℓArray[ℓ,:,:])
        vecpCℓMCMC = zeros(floor(Int,inumber*0.5*(inumber+1)))
        LinearAlgebra.mul!(vecpCℓMCMC, L, vecCℓMCMC)
    end

    for ℓidx in eachindex(Cl)
        vecpCℓMCMC[ℓidx,:,:] ~ MvNormal(vecpCℓData[ℓidx,:,:], Cov.CℓCovariance[ℓidx,:,:])
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
@everywhere Covaₗₘ = InstantiateEvaluateCovariance(CentralCℓ, ConvolvedDensity, CosmologicalGrid, "Lensing",
"Lensing")
@everywhere CovCℓ = InstantiateEvaluateCovariance(Covaₗₘ)

@everywhere model = gdemo(CℓData, CovCℓ, ConvolvedDensity, EuclidBias, EuclidIA,
CosmologicalGrid)

chains = sample(model, MH(), MCMCDistributed(), 100, n_process; save_state = true)
write("Turing/chains/first_chain-file.jls", chains)