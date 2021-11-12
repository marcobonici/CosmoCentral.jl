```@setup tutorial
using Plots; gr()
Plots.reset_defaults()
using CosmoCentral
using JSON
using LaTeXStrings
using CairoMakie
using Makie
using LinearAlgebra
using FisherPlot
run(`wget https://zenodo.org/record/5270335/files/forecast_pmm.tar.xz\?download=1`);
run(`mv forecast_pmm.tar.xz\?download\=1 forecast_pmm.tar.xz`);
run(`tar xvf forecast_pmm.tar.xz`);
run(`rm -rf forecast_pmm.tar.xz`);
MultipolesArrayTemp = CosmoCentral.LogSpaced(10.,3000., 101)
MultipolesArray = zeros(100)
MultipolesWidths = CosmoCentral.Difference(MultipolesArrayTemp)
for i in 1:100
    MultipolesArray[i] = (MultipolesArrayTemp[i+1]+MultipolesArrayTemp[i])/2
end

steps = Array([0.00625, 0.01250, 0.01875, 0.02500, 0.03750, 0.05000, 0.10000])
cosmogrid = CosmoCentral.CosmologicalGrid(ZArray =
Array(LinRange(0.001, 4., 500)), KArray = CosmoCentral.LogSpaced(1e-5, 50., 1000),
ℓBinCenters = MultipolesArray, ℓBinWidths = MultipolesWidths)
ProbesDict = JSON.parsefile(pwd()*"/../../input_files/AngularNew.json")
CosmoDict = JSON.parsefile(pwd()*"/../../input_files/Cosmology.json")
ForecastContainer = CosmoCentral.InitializeForecastContainer(CosmoDict, ProbesDict,
cosmogrid, steps)
CosmoCentral.CreateDirectoriesForecast!(ForecastContainer, pwd()*"/test_forecast/")
PathInputPmm = pwd()*"/forecast_pmm/PowerSpectrum/"
PathOutputCℓ = pwd()*"/test_forecast/Angular/"
PathOutput = pwd()*"/test_forecast"
PathCentralCℓ = pwd()*"/test_forecast/Angular/dvar_central_step_0/cl"
Path∂Cℓ = pwd()*"/test_forecast/Derivative"
CosmoCentral.ForecastCℓ!(ForecastContainer, cosmogrid, PathInputPmm, PathOutputCℓ)
CosmoCentral.Forecast∂Cℓ!(ForecastContainer, PathOutput, PathOutputCℓ, steps)
FisherWL = CosmoCentral.ForecastFisherαβ(ForecastContainer ,PathCentralCℓ, Path∂Cℓ,
cosmogrid, "Lensing", "Test")
FisherGC = CosmoCentral.ForecastFisherαβ(ForecastContainer ,PathCentralCℓ, Path∂Cℓ,
cosmogrid, "PhotometricClustering", "Test")
FisherXC_WL_GC = CosmoCentral.ForecastFisherαβ(ForecastContainer ,PathCentralCℓ, Path∂Cℓ,
cosmogrid, "Lensing", "PhotometricClustering", "Test")
FisherGCplusWL = CosmoCentral.SumFisher(FisherWL, FisherGC)
FisherGCNew = CosmoCentral.RearrangeFisherℓ(FisherGC, cosmogrid, 10, 3000)
CosmoCentral.SelectMatrixAndMarginalize!(FisherGCNew.ParametersList, FisherGCNew)
FisherGCNewplusWL = CosmoCentral.SumFisher(FisherWL, FisherGCNew)
FisherXC_WL_GCNew = CosmoCentral.RearrangeFisherℓ(FisherXC_WL_GC, cosmogrid, 10, 3000)
FisherWL_Add = CosmoCentral.RearrangeFisherℓ(FisherWL, cosmogrid, 3000, 5000)
FisherFinal = CosmoCentral.SumFisher(FisherWL_Add, FisherXC_WL_GCNew)

BigLaTeXArray = [L"\Omega_\mathrm{M}", L"\Omega_\mathrm{B}", L"H_0", L"n_s", L"\sigma_8", L"M_\nu", L"w_a", L"w_0"]
pars_list = ["ΩM", "ΩB", "H0", "ns", "σ8", "Mν", "wa", "w0"]
central_values =[0.32, 0.05, 67., 0.96, 0.816, 0.06, 0., -1.0]

marg_errs = []
for par in pars_list
    append!(marg_errs, FisherWL.MarginalizedErrors[par])
end

noto_sans = "Computer Modern"#assetpath("/usr/share/fonts/computer_modern", "NewCM10-Regular.otf")

dimticklabel = 50
sidesquare = 400
ylabelsize = 39

PlotPars = Dict("sidesquare" => sidesquare,
"dimticklabel" => dimticklabel,
"parslabelsize" => 80,
"textsize" => 80,
"PPmaxlabelsize" => 60,
"font" => noto_sans,
"xticklabelrotation" => 45.)

limits = zeros(8,2)
ticks = zeros(8,2)
for i in 1:8
    limits[i,1] = central_values[i]-3marg_errs[i]*1.3
    limits[i,2] = central_values[i]+3marg_errs[i]*1.3
    ticks[i,1] = central_values[i]-3marg_errs[i]
    ticks[i,2] = central_values[i]+3marg_errs[i]
end

probes = [L"\mathrm{WL}", L"\mathrm{GC}_\mathrm{ph}",
L"\mathrm{WL}\,+\,\mathrm{GC}_\mathrm{ph}",
L"\mathrm{WL}\,+\,\mathrm{GC}_\mathrm{ph}\,+\,\mathrm{XC}"]
colors = ["deepskyblue3", "darkorange1", "green", "red"]
```

# Fisher Forecast

Fisher Matrices can be used to forecast parameter uncertainties with few computational
resources. They are encapsulated in the following struct
```@docs
CosmoCentral.Fisherαβ
```
We implement two different approaches to evaluate the Fisher Matrix:

- A field approach, where we consider the ``a_{\ell m}`` to be the observable
- An estimator approach, where we consider the ``C_\ell``'s to be the observable

If correctly implemented, the two approaches give the same result ([Hamimeche & Lewis 2008](https://arxiv.org/abs/0801.0554), [Carron 2012](https://arxiv.org/abs/1204.4724)).

## Field approach

The Fisher Matrix in the Field approach is evaluated by the following method

```@docs
CosmoCentral.ForecastFisherαβ(PathCentralCℓ::String, Path∂Cℓ::String,
InputList::Vector{Dict{String, Vector{Any}}},
CosmologicalGrid::CosmoCentral.CosmologicalGrid)
```

## Estimator approach

The Fisher Matrix in the Estimator approach is evaluated by the following method

```@docs
CosmoCentral.ForecastFisherαβ(PathCentralCℓ::String, Path∂Cℓ::String,
InputList::Vector{Dict{String, Vector{Any}}},
CosmologicalGrid::CosmoCentral.CosmologicalGrid, ciccio::String)
```

## Plotting Fisher Matrix
Here we show plots for some Fisher Matrices (this code will probably soon released in a
separate package).
```@example tutorial
canvas = FisherPlot.PrepareCanvas(BigLaTeXArray, central_values, limits, ticks, probes, colors,
PlotPars)
FisherPlot.PaintCorrMattrix!(canvas, central_values,
CosmoCentral.SelectCorrelationMatrix(FisherWL, pars_list), "deepskyblue3")
FisherPlot.PaintCorrMattrix!(canvas, central_values,
CosmoCentral.SelectCorrelationMatrix(FisherGCNew, pars_list), "darkorange1")
FisherPlot.PaintCorrMattrix!(canvas, central_values,
CosmoCentral.SelectCorrelationMatrix(FisherGCNewplusWL, pars_list), "green")
FisherPlot.PaintCorrMattrix!(canvas, central_values,
CosmoCentral.SelectCorrelationMatrix(FisherFinal, pars_list), "red")
canvas
```