# Load Turing and MCMCChains.
using Turing
using CosmoCentral
using Distributed
using Dates
using ClusterManagers
using LinearAlgebra

n_process = 100
n_step = 400
#add  n_processes on the long queue
for i in 1:n_process
    time = string(Dates.now())
    ClusterManagers.addprocs_lsf(1; bsub_flags = `-q long -o Turing/outfiles/$time.out -e Turing/outfiles/$time.err -M 3000 -R 'hname!=teo26.ge.infn.it && hname!=teo25.ge.infn.it && hname!=teo27.ge.infn.it && hname!=teo28.ge.infn.it && hname!=g2farm7.ge.infn.it && hname!=totem04.ge.infn.it && hname!=totem07.ge.infn.it && hname!=totem08.ge.infn.it &&  hname!=g2farm8.ge.infn.it && hname!=teo22.ge.infn.it && hname!=teo24.ge.infn.it && hname!=teo28.ge.infn.it'`)
end

@everywhere using Distributed, Turing, LinearAlgebra
@everywhere using CosmoCentral

@everywhere w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()
println("Loaded central cosmology")

@everywhere @model function gdemo(vecpCℓData, Cov, ConvolvedDensity, EuclidBias, EuclidIA,
    CosmologicalGrid)
    w₀ ~ Uniform(-1.32, -0.72)
    wₐ ~ Uniform(-0.68, 0.72)
    ns ~ Uniform(0.92, 1.0)
    ΩM ~ Uniform(0.29, 0.35)
    ΩB ~ Uniform(0.04, 0.06)  
    σ8 ~ Uniform(0.6,  1.0)
    H0 ~ Uniform(61.,  73)

    w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()
    w0waCDMCosmology.w0 = w₀
    w0waCDMCosmology.wa = wₐ
    w0waCDMCosmology.ns = ns
    w0waCDMCosmology.ΩM = ΩM
    w0waCDMCosmology.ΩB = ΩB
    w0waCDMCosmology.σ8 = σ8
    w0waCDMCosmology.H0 = H0
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

chains = sample(model, MH(), MCMCDistributed(), n_step, n_process; save_state = true)
write("chain_one.jls", chains)
