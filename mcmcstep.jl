using Revise
using CosmoCentral
w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()

MultipolesArrayTemp = CosmoCentral.LogSpaced(10.,3000., 101)
MultipolesArray = zeros(100)
#MultipolesWidths = vcat(CosmoCentral.Difference(MultipolesArrayTemp), ones(2000))
MultipolesWidths = CosmoCentral.Difference(MultipolesArrayTemp)
for i in 1:100
    MultipolesArray[i] = (MultipolesArrayTemp[i+1]+MultipolesArrayTemp[i])/2
end
CosmologicalGrid = CosmoCentral.CosmologicalGrid(ZArray = LinRange(0.001, 4., 500), KArray = CosmoCentral.LogSpaced(1e-5, 50., 1000), ℓBinCenters = MultipolesArray)

AnalitycalDensity = CosmoCentral.AnalitycalDensity()
CosmoCentral.NormalizeAnalitycalDensity!(AnalitycalDensity)

InstrumentResponse = CosmoCentral.InstrumentResponse()

ConvolvedDensity = CosmoCentral.ConvolvedDensity(DensityGridArray = ones(10, length(CosmologicalGrid.ZArray)))
CosmoCentral.NormalizeConvolvedDensity!(ConvolvedDensity, AnalitycalDensity, InstrumentResponse, CosmologicalGrid)

CosmoCentral.ComputeConvolvedDensityGrid!(CosmologicalGrid, ConvolvedDensity, AnalitycalDensity, InstrumentResponse)

EuclidBias = CosmoCentral.EuclidBias()
EuclidIA = CosmoCentral.ExtendedNLIA()

CℓMCMC = CosmoCentral.EvaluateCℓMCMCStep(w0waCDMCosmology, ConvolvedDensity, EuclidBias, EuclidIA, CosmologicalGrid);
