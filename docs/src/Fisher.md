```@setup tutorial
using Plots; gr()
Plots.reset_defaults()
using CosmoCentral
using JSON
using LaTeXStrings
using CairoMakie
using Makie
using LinearAlgebra
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
CosmologicalGrid = CosmoCentral.CosmologicalGrid(ZArray =
Array(LinRange(0.001, 4., 500)), KArray = CosmoCentral.LogSpaced(1e-5, 50., 1000),
ℓBinCenters = MultipolesArray, ℓBinWidths = MultipolesWidths)
ProbesDict = JSON.parsefile(pwd()*"/../../input_files/AngularNew.json")
CosmoDict = JSON.parsefile(pwd()*"/../../input_files/Cosmology.json")
ForecastContainer = CosmoCentral.InitializeForecastContainer(CosmoDict, ProbesDict,
CosmologicalGrid, steps)
CosmoCentral.CreateDirectoriesForecast!(ForecastContainer, pwd()*"/test_forecast/")
PathInputPmm = pwd()*"/forecast_pmm/PowerSpectrum/"
PathOutputCℓ = pwd()*"/test_forecast/Angular/"
PathOutput = pwd()*"/test_forecast"
PathCentralCℓ = pwd()*"/test_forecast/Angular/dvar_central_step_0/cl"
Path∂Cℓ = pwd()*"/test_forecast/Derivative"
CosmoCentral.ForecastCℓ!(ForecastContainer, CosmologicalGrid, PathInputPmm, PathOutputCℓ)
CosmoCentral.Forecast∂Cℓ!(ForecastContainer, PathOutput, PathOutputCℓ, steps)
FisherWL = CosmoCentral.ForecastFisherαβ(ForecastContainer ,PathCentralCℓ, Path∂Cℓ,
CosmologicalGrid, "Lensing", "Test")
FisherGC = CosmoCentral.ForecastFisherαβ(ForecastContainer ,PathCentralCℓ, Path∂Cℓ,
CosmologicalGrid, "PhotometricClustering", "Test")
FisherXC_WL_GC = CosmoCentral.ForecastFisherαβ(ForecastContainer ,PathCentralCℓ, Path∂Cℓ,
CosmologicalGrid, "Lensing", "PhotometricClustering", "Test")
FisherGCplusWL = CosmoCentral.SumFisher(FisherWL, FisherGC)
FisherGCNew = CosmoCentral.RearrangeFisherℓ(FisherGC, CosmologicalGrid, 10, 3000)
CosmoCentral.SelectMatrixAndMarginalize!(FisherGCNew.ParametersList, FisherGCNew)
FisherGCNewplusWL = CosmoCentral.SumFisher(FisherWL, FisherGCNew)
FisherXC_WL_GCNew = CosmoCentral.RearrangeFisherℓ(FisherXC_WL_GC, cosmogrid, 10, 3000)
FisherWL_Add = CosmoCentral.RearrangeFisherℓ(FisherWL, cosmogrid, 3000, 5000)
FisherFinal = CosmoCentral.SumFisher(FisherWL_Add, FisherXC_WL_GCNew)

function EllipseParametrization(a::Float64, b::Float64, θ::Float64)
    t = LinRange(0,2π, 200)
    x = Array(a .* sin.(θ) .* cos.(t) + b .* cos.(θ) .* sin.(t))
    y = Array(a .* cos.(θ) .* cos.(t) - b .* sin.(θ) .* sin.(t))
    return x, y
end

function Gaussian(μ::Float64, σ::Float64, x)
    return 1/(sqrt(2π*σ^2))*exp(-0.5*(x-μ)^2/σ^2)
end

function EllipseParameters(covmatrix::Matrix{Float64}, i::Int64, j::Int64)
    σi = sqrt(covmatrix[i,i])
    σj = sqrt(covmatrix[j,j])
    σij = covmatrix[i,j]
    θ = (atan(2σij,(σi^2-σj^2)))/2
    a = sqrt((σi^2+σj^2)/2+sqrt(((σi^2-σj^2)^2)/4+σij^2))
    if i==j
        b = 0.
    else
        b = sqrt((σi^2+σj^2)/2-sqrt(((σi^2-σj^2)^2)/4+σij^2))
    end
    return σi, σj, a, b, θ
end

marg_errs = []
for (key, value) in FisherGCNew.MarginalizedErrors
    append!(marg_errs, value)
end

BigLaTeXArray = [L"\Omega_\mathrm{B}", L"M_\nu", L"w_a", L"n_s", L"\Omega_\mathrm{M}", L"\sigma_8", L"w_0", L"H_0"]
pars_list = ["ΩB", "Mν", "wa", "ns",  "ΩM", "σ8", "w0", "H0"]
central_values =[0.05, 0.06, 0., 0.96, 0.32, 0.816, -1.0, 67.]

BigLaTeXArray = [L"\Omega_\mathrm{M}", L"\Omega_\mathrm{B}", L"H_0", L"n_s", L"\sigma_8", L"M_\nu", L"w_a", L"w_0"]
pars_list = ["ΩM", "ΩB", "H0", "ns", "σ8", "Mν", "wa", "w0"]
central_values =[0.32, 0.05, 67., 0.96, 0.816, 0.06, 0., -1.0]

noto_sans = assetpath("/usr/share/fonts/computer_modern", "NewCM10-Regular.otf")

dimticklabel = 50#36
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

function PrepareCanvas(LaTeX_array, central_values, limits, ticks, probes, colors, PlotPars)
    matrix_dim = length(LaTeX_array)
    #TODO: add the textsize to PlotPars
    figure = Makie.Figure(textsize = 40, font = PlotPars["font"])
    for i in 1:matrix_dim
        for j in 1:i
            if i == j
                ax = Axis(figure[i,i],
                    width = PlotPars["sidesquare"], height = PlotPars["sidesquare"],
                    xticklabelsize = PlotPars["dimticklabel"],
                    yticklabelsize = PlotPars["dimticklabel"], yaxisposition = (:right),
                    xlabel = L"%$((LaTeX_array)[i])", xlabelsize = PlotPars["parslabelsize"],
                    ylabel = L"P/P_\mathrm{max}",ylabelsize = PlotPars["PPmaxlabelsize"], yticks = [0,1],
                    xticklabelrotation = PlotPars["xticklabelrotation"],
                    xticks = ([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]],
                        [string(myi) for myi in round.([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]], sigdigits = 3)]))
                Makie.ylims!(ax, (-0.0,1.05))
                Makie.xlims!(ax, (limits[i,1],limits[i,2]))
                Makie.hideydecorations!(ax, ticks = false, ticklabels = false, label = false)
                if i != matrix_dim
                    ax.alignmode = Mixed(right = MakieLayout.Protrusion(0), bottom = MakieLayout.Protrusion(0), top= MakieLayout.Protrusion(0))
                    hidexdecorations!(ax, ticks = true, ticklabels = true,  label = true)
                else
                    hidexdecorations!(ax, ticks = false, ticklabels = false,  label = false)
                end
            else
                ax = Axis(figure[i,j], width = PlotPars["sidesquare"], height = PlotPars["sidesquare"],
                    xticklabelsize = PlotPars["dimticklabel"], yticklabelsize = PlotPars["dimticklabel"],
                    ylabel = L"%$(LaTeX_array[i])", xlabel = L"%$((LaTeX_array)[j])",
                    ylabelsize = PlotPars["parslabelsize"], xlabelsize = PlotPars["parslabelsize"], xticklabelrotation = PlotPars["xticklabelrotation"],
                    yticks = ([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]],
                        [string(myi) for myi in round.([ticks[i,1], 0.5*(ticks[i,1]+ticks[i,2]), ticks[i,2]], sigdigits = 3)]),
                    xticks = ([ticks[j,1], 0.5*(ticks[j,1]+ticks[j,2]), ticks[j,2]],
                        [string(myi) for myi in round.([ticks[j,1], 0.5*(ticks[j,1]+ticks[j,2]), ticks[j,2]], sigdigits = 3)]),
                yticklabelpad=8)
                Makie.ylims!(ax, (limits[i,1],limits[i,2]))
                Makie.xlims!(ax, (limits[j,1],limits[j,2]))
                if i == matrix_dim
                    hidexdecorations!(ax, ticks = false, ticklabels = false,  label = false)
                else
                    hidexdecorations!(ax, ticks = true, ticklabels = true,  label = true)
                end
                if j == 1
                    hideydecorations!(ax, ticks = false, ticklabels = false,  label = false)
                    Legend(figure[1,matrix_dim],
                    [PolyElement(color = color, strokecolor = color, strokewidth = 1) for color in colors],
                    probes,
                    tellheight = false, tellwidth = false, rowgap = 10,
                    halign = :right, valign = :top, framecolor = :black, labelsize =55, patchsize = (70, 40), framevisible = true)
                else
                    hideydecorations!(ax, ticks = true, ticklabels = true,  label = true)
                    ax.alignmode = Mixed(right = MakieLayout.Protrusion(0), bottom = MakieLayout.Protrusion(0), top = MakieLayout.Protrusion(0))
                end
            end
        end
    end
    Makie.resize!(figure.scene, figure.layout.layoutobservables.reportedsize[]...)
    return figure
end

probes = [L"\mathrm{WL}", L"\mathrm{GC}_\mathrm{ph}", L"\mathrm{WL}\,+\,\mathrm{GC}_\mathrm{ph}", L"\mathrm{WL}\,+\,\mathrm{GC}_\mathrm{ph}\,+\,\mathrm{XC}"]
colors = ["deepskyblue3", "darkorange1", "green", "red"]

function DrawGaussian!(canvas, σ, i, central, color)
    ax = canvas[i,i]
    x = Array(LinRange(-4σ+central,4σ+central, 200))
    Makie.lines!(ax, x, Gaussian.(central, σ, x)./Gaussian.(central, σ, central), color = color, linewidth = 4)
    Makie.band!(ax, x, 0, Gaussian.(central, σ, x)./Gaussian.(central, σ, central) , color=(color, 0.2))
    x = Array(LinRange(-σ+central,σ+central, 200))
    Makie.band!(ax, x, 0, Gaussian.(central, σ, x)./Gaussian.(central, σ, central) , color=(color, 0.4))
end 

function DrawEllipse!(canvas, i, j, x, y, central_values, color)
    ax = canvas[i,j]
    
    Makie.lines!(ax, x .+ central_values[j], y .+ central_values[i], color = color, linewidth = 4)
    Makie.lines!(ax, 3x .+ central_values[j], 3y .+ central_values[i], color = color, linewidth = 4)

    Makie.band!(ax, x .+ central_values[j], 0, y .+ central_values[i], color=(color, 0.4))
    Makie.band!(ax, 3x .+ central_values[j], 0, 3y .+ central_values[i], color=(color, 0.2))
end

function PaintCorrMattrix!(canvas, central_values, corr_matrix, color)
    for i in 1:length(central_values)
        for j in 1:i
            if i == j
                DrawGaussian!(canvas, sqrt(corr_matrix[i,i]), i, central_values[i], color)
            else
                σi, σj, a, b, θ = EllipseParameters(corr_matrix, i,j)
                x,y = EllipseParametrization(a, b, θ)
                DrawEllipse!(canvas, i, j, x, y, central_values, color)
            end
        end
    end
end


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
canvas = PrepareCanvas(BigLaTeXArray, central_values, limits, ticks, probes, colors,
PlotPars)
PaintCorrMattrix!(canvas, central_values,
CosmoCentral.SelectCorrelationMatrix(FisherWL, pars_list), "deepskyblue3")
PaintCorrMattrix!(canvas, central_values,
CosmoCentral.SelectCorrelationMatrix(FisherGCNew, pars_list), "darkorange1")
PaintCorrMattrix!(canvas, central_values,
CosmoCentral.SelectCorrelationMatrix(FisherGCNewplusWL, pars_list), "green")
PaintCorrMattrix!(canvas, central_values,
CosmoCentral.SelectCorrelationMatrix(FisherFinal, pars_list), "red")
canvas
```