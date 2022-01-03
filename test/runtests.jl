using CosmoCentral
using Test
using QuadGK
using NumericalIntegration
using PyCall
using Interpolations
using JSON

numpy = pyimport("numpy")

w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()
AnalitycalDensity = CosmoCentral.AnalitycalDensity()
InstrumentResponse = CosmoCentral.InstrumentResponse()
ConvolvedDensity = CosmoCentral.ConvolvedDensity()
CosmologicalGrid  = CosmoCentral.CosmologicalGrid(
ZArray=Array(LinRange(0.001, 4.0, 500)))
BackgroundQuantities = CosmoCentral.BackgroundQuantities(HZArray=
zeros(length(CosmologicalGrid.ZArray)),
χZArray=zeros(length(CosmologicalGrid.ZArray)))
GCWeightFunction = CosmoCentral.GCWeightFunction(WeightFunctionArray =
zeros(length(ConvolvedDensity.DensityNormalizationArray),
length(CosmologicalGrid.ZArray)))
WLWeightFunction = CosmoCentral.WLWeightFunction(WeightFunctionArray =
zeros(length(ConvolvedDensity.DensityNormalizationArray),
length(CosmologicalGrid.ZArray)),
LensingEfficiencyArray = zeros(length(
ConvolvedDensity.DensityNormalizationArray),
length(CosmologicalGrid.ZArray)))
AngularCoefficients = CosmoCentral.Cℓ(
CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters),
length(GCWeightFunction.WeightFunctionArray[:, 1]),
length(GCWeightFunction.WeightFunctionArray[:, 1])))
classyParams = CosmoCentral.Initializeclassy(w0waCDMCosmology)
input_path_pmm = pwd()*"/forecast_pmm/PowerSpectrum/dvar_central_step_0/p_mm"
input_path_Cℓ = pwd()*"/forecast_pmm/Angular/dvar_central_step_0/cl"
input_path_Forecast = pwd()
CosmoCentral.classy.Class()


run(`wget https://zenodo.org/record/5270335/files/forecast_pmm.tar.xz\?download=1`);
run(`mv forecast_pmm.tar.xz\?download\=1 forecast_pmm.tar.xz`);
run(`tar xvf forecast_pmm.tar.xz`);

@testset "Evaluation of background quantities" begin
    test_E_z = CosmoCentral.ComputeAdimensionalHubbleFactor(0., w0waCDMCosmology)
    @test test_E_z == 1.
    test_H_z = CosmoCentral.ComputeHubbleFactor(0., w0waCDMCosmology)
    @test test_H_z == w0waCDMCosmology.H0
    test_r_z = CosmoCentral.Computeχ(0., w0waCDMCosmology)
    @test test_r_z == 0.
    test_r_array = zeros(length(CosmologicalGrid.ZArray))
    test_H_array = zeros(length(CosmologicalGrid.ZArray))
    for (idxz, zvalue) in enumerate(CosmologicalGrid.ZArray)
        test_r_array[idxz] = CosmoCentral.Computeχ(zvalue,
        w0waCDMCosmology)
        test_H_array[idxz] = CosmoCentral.ComputeHubbleFactor(zvalue,
        w0waCDMCosmology)
    end
    CosmoCentral.ComputeBackgroundQuantitiesGrid!(CosmologicalGrid,
    BackgroundQuantities, w0waCDMCosmology)
    @test test_r_array ==BackgroundQuantities.χZArray
    @test test_H_array ==BackgroundQuantities.HZArray
end

@testset "Check the normalization of density function" begin
    CosmoCentral.NormalizeAnalitycalDensity!(AnalitycalDensity)
    int, err = QuadGK.quadgk(x -> CosmoCentral.ComputeDensity(x, AnalitycalDensity),
    AnalitycalDensity.ZMin, AnalitycalDensity.ZMax, rtol=1e-12)
    @test isapprox(int, AnalitycalDensity.SurfaceDensity, atol=1e-9)
end

@testset "Check the normalization of convolved density function" begin
    test_normalization = zeros(length(ConvolvedDensity.ZBinArray)-1)
    CosmoCentral.NormalizeConvolvedDensity!(ConvolvedDensity, AnalitycalDensity,
    InstrumentResponse, CosmologicalGrid)
    for idx in 1:length(test_normalization)
        int, err = QuadGK.quadgk(x ->
        CosmoCentral.ComputeConvolvedDensity(x, idx,
        ConvolvedDensity, AnalitycalDensity, InstrumentResponse),
        AnalitycalDensity.ZMin,
        AnalitycalDensity.ZMax, rtol=1e-12)
        test_normalization[idx] = int
    end
    @test isapprox(test_normalization, ones(length(test_normalization)), atol=1e-12)
end

@testset "Check the computation of the convolved density function on grid" begin
    test_array = zeros(Float64,length(ConvolvedDensity.ZBinArray)-1,
    length(CosmologicalGrid.ZArray))
    CosmoCentral.ComputeConvolvedDensityGrid!(CosmologicalGrid, ConvolvedDensity,
    AnalitycalDensity, InstrumentResponse)
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            test_array[idx_ZBinArray, idx_ZArray] =
            CosmoCentral.ComputeConvolvedDensity(
            CosmologicalGrid.ZArray[idx_ZArray],
            idx_ZBinArray,
            ConvolvedDensity,
            AnalitycalDensity,
            InstrumentResponse)
        end
    end
    @test isapprox(test_array, ConvolvedDensity.DensityGridArray, atol=1e-12)
end

@testset "Check the functions defined in MathUtils" begin
    minarray = 1e-5
    maxarray = 10.
    n = 100
    benchmark_numpy = numpy.logspace(log10(minarray), log10(maxarray), n)
    test_logspace = CosmoCentral.LogSpaced(minarray, maxarray, n)
    @test isapprox(benchmark_numpy, test_logspace, atol=1e-12)
    array = [1.,2.,3.]
    @test 1 == CosmoCentral.BinSearch(1.5, array)
    x = Array(LinRange(0., 10., 100))
    y = 11.2 .*x .+0.4
    c, m = CosmoCentral.CustomRegression(x, y)
    @test isapprox(c, 0.4, atol=1e-12)
    @test isapprox(m, 11.2, atol=1e-12)
end

@testset "Check the Piecewise bias evaluation" begin
    z = 0.4
    test_array = zeros(length(CosmologicalGrid.ZArray))
    bias = CosmoCentral.ComputeBias(z, GCWeightFunction.BiasKind,
    ConvolvedDensity)
    @test isapprox(1.0997727037892875, bias, atol=1e-12)
    for (idxz, zvalue) in enumerate(CosmologicalGrid.ZArray)
        test_array[idxz] = CosmoCentral.ComputeBias(zvalue,
        GCWeightFunction.BiasKind, ConvolvedDensity)
    end
    CosmoCentral.ComputeBiasGrid!(CosmologicalGrid,
        GCWeightFunction, ConvolvedDensity)
    @test isapprox(test_array, GCWeightFunction.BiasArray[1,:], atol=1e-10)
end

@testset "Check the Weight function evaluation" begin
    test_gc = zeros(length(ConvolvedDensity.ZBinArray)-1,
    length(CosmologicalGrid.ZArray))
    test_le = zeros(length(ConvolvedDensity.ZBinArray)-1,
    length(CosmologicalGrid.ZArray))
    test_wl = zeros(length(ConvolvedDensity.ZBinArray)-1,
    length(CosmologicalGrid.ZArray))
    for (idxz, zvalue) in enumerate(CosmologicalGrid.ZArray)
        for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
            test_gc[idx_ZBinArray, idxz] =
            CosmoCentral.ComputeWeightFunction(zvalue, idx_ZBinArray,
            ConvolvedDensity, AnalitycalDensity, InstrumentResponse,
            w0waCDMCosmology, GCWeightFunction)
            test_le[idx_ZBinArray, idxz] =
            CosmoCentral.ComputeLensingEfficiency(zvalue, idx_ZBinArray,
            ConvolvedDensity, AnalitycalDensity, InstrumentResponse,
            w0waCDMCosmology, CosmologicalGrid, WLWeightFunction)
            test_wl[idx_ZBinArray, idxz] =
            CosmoCentral.ComputeWeightFunction(zvalue, idx_ZBinArray,
                ConvolvedDensity, AnalitycalDensity,
                InstrumentResponse, w0waCDMCosmology, CosmologicalGrid,
                WLWeightFunction)
        end
    end
    CosmoCentral.ComputeWeightFunctionGrid!(GCWeightFunction,
        ConvolvedDensity,
        CosmologicalGrid,
        BackgroundQuantities,
        w0waCDMCosmology)
    CosmoCentral.ComputeLensingEfficiencyGrid!(
    WLWeightFunction, AnalitycalDensity,
    InstrumentResponse, ConvolvedDensity,
    CosmologicalGrid,
    BackgroundQuantities,
    w0waCDMCosmology,
    CosmoCentral.StandardLensingEfficiency())
    CosmoCentral.ComputeWeightFunctionGrid!(WLWeightFunction,
        ConvolvedDensity,
        CosmologicalGrid,
        BackgroundQuantities,
        w0waCDMCosmology)
    @test isapprox(test_gc, GCWeightFunction.WeightFunctionArray, rtol=1e-9)
    @test isapprox(test_wl, WLWeightFunction.WeightFunctionArray, rtol=1e-9)
    @test isapprox(test_le, WLWeightFunction.LensingEfficiencyArray, rtol=1e-9)
    CosmoCentral.ComputeIntrinsicAlignmentGrid!(CosmologicalGrid, WLWeightFunction, ConvolvedDensity, BackgroundQuantities, w0waCDMCosmology)
    CosmoCentral.ComputeWeightFunctionGrid!(WLWeightFunction, ConvolvedDensity, CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
end


@testset "Check the Power Spectrum evaluated over the Limber Grid" begin
    ℓBinCenters = Array(LinRange(10.5, 2999.5, 2990))
    PowerSpectrum, BackgroundQuantitiesLoaded, CosmologicalGrid =
    CosmoCentral.ReadPowerSpectrumBackground(input_path_pmm, ℓBinCenters, Array(ones(2990)))
    CosmoCentral.InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
    BackgroundQuantitiesLoaded, PowerSpectrum, CosmoCentral.BSplineCubic())
    CosmoCentral.ComputeLimberArray!(CosmologicalGrid,
    BackgroundQuantitiesLoaded)
    test_k_limber = (CosmologicalGrid.ℓBinCenters[1]+0.5) /
    BackgroundQuantitiesLoaded.χZArray[1]
    @test test_k_limber == CosmologicalGrid.KLimberArray[1, 1]
    test_Omega_cdm = w0waCDMCosmology.ΩM-w0waCDMCosmology.ΩB-
    w0waCDMCosmology.Mν/(93.14*(w0waCDMCosmology.H0/100)^2)
    @test test_Omega_cdm == classyParams.classyParamsDict["Omega_cdm"]
    x = Array(LinRange(log10(first(CosmologicalGrid.KArray)),
    log10(last(CosmologicalGrid.KArray)), length(CosmologicalGrid.KArray)))
    y = Array(LinRange(first(CosmologicalGrid.ZArray), last(CosmologicalGrid.ZArray),
    length(CosmologicalGrid.ZArray)))
    InterpPmm = Interpolations.interpolate(
    log10.(PowerSpectrum.PowerSpectrumNonlinArray),
    BSpline(Cubic(Line(OnGrid()))))
    InterpPmm = Interpolations.extrapolate(InterpPmm, Line())
    test_power_spectrum = 10 .^InterpPmm(log10(
    CosmologicalGrid.KLimberArray[1, 1]), CosmologicalGrid.ZArray[1])
    CosmoCentral.InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
    BackgroundQuantitiesLoaded, PowerSpectrum, CosmoCentral.BSplineCubic())
    @test isapprox(test_power_spectrum,
    PowerSpectrum.InterpolatedPowerSpectrum[1, 1], rtol=1e-2)
    CosmoCentral.WritePowerSpectrumBackground(PowerSpectrum,
    BackgroundQuantitiesLoaded, CosmologicalGrid, "new_p_mm")
    NewPowerSpectrum, NewBackgroundQuantitiesLoaded, NewCosmologicalGrid =
    CosmoCentral.ReadPowerSpectrumBackground("new_p_mm", ℓBinCenters, Array(ones(2990)))
    @test isapprox(PowerSpectrum.PowerSpectrumNonlinArray,
    NewPowerSpectrum.PowerSpectrumNonlinArray, rtol=1e-6)
    @test isapprox(PowerSpectrum.PowerSpectrumLinArray,
    NewPowerSpectrum.PowerSpectrumLinArray, rtol=1e-6)
    @test isapprox(BackgroundQuantitiesLoaded.HZArray,
    NewBackgroundQuantitiesLoaded.HZArray, rtol=1e-6)
    @test isapprox(BackgroundQuantitiesLoaded.χZArray,
    NewBackgroundQuantitiesLoaded.χZArray, rtol=1e-6)
    PowerSpectrumDierckx, BackgroundQuantitiesLoaded, CosmologicalGrid =
    CosmoCentral.ReadPowerSpectrumBackground(input_path_pmm, ℓBinCenters, Array(ones(2990)))
    CosmoCentral.InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
    BackgroundQuantitiesLoaded, PowerSpectrumDierckx,
    CosmoCentral.RectBivSplineDierckx())
    @test isapprox(PowerSpectrumDierckx.InterpolatedPowerSpectrum[1, 1],
    PowerSpectrum.InterpolatedPowerSpectrum[1, 1], rtol=1e-2)
end

"""
@testset "Test Angular coefficients evaluation" begin
    CℓLoaded = CosmoCentral.ReadCℓ(input_path_Cℓ, "PhotometricGalaxy_PhotometricGalaxy")
    MultipolesArrayTemp = CosmoCentral.LogSpaced(10.,3000., 101)
    MultipolesArray = zeros(100)
    MultipolesWidths = CosmoCentral.Difference(MultipolesArrayTemp)
    for i in 1:100
        MultipolesArray[i] = (MultipolesArrayTemp[i+1]+MultipolesArrayTemp[i])/2
    end
    ℓBinCenters = MultipolesArray
    PowerSpectrum, BackgroundQuantities, CosmologicalGrid =
    CosmoCentral.ReadPowerSpectrumBackground(input_path_pmm, ℓBinCenters, Array(ones(2990)))
    CosmoCentral.ComputeLimberArray!(CosmologicalGrid, BackgroundQuantities)
    CosmoCentral.InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
    BackgroundQuantities, PowerSpectrum, CosmoCentral.BSplineCubic())
    Cℓ = CosmoCentral.Cℓ(CℓArray = zeros(length(CosmologicalGrid.ℓBinCenters),
    length(GCWeightFunction.WeightFunctionArray[:, 1]),
    length(GCWeightFunction.WeightFunctionArray[:, 1])))
    CosmoCentral.ComputeCℓ!(Cℓ, GCWeightFunction, GCWeightFunction, BackgroundQuantities,
    w0waCDMCosmology, CosmologicalGrid, PowerSpectrum, CosmoCentral.CustomSimpson())
    @test isapprox(CℓLoaded.CℓArray,
    Cℓ.CℓArray, rtol=1e-6)
    CosmoCentral.WriteCℓ!("PhotometricGalaxy_PhotometricGalaxy", Cℓ, "new_cl")
    CℓReloaded = CosmoCentral.ReadCℓ("new_cl", "PhotometricGalaxy_PhotometricGalaxy")
    @test isapprox(CℓReloaded.CℓArray, Cℓ.CℓArray, rtol=1e-9)
end
"""

@testset "Test Fisher Forecast: Integration Test" begin
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
    ProbesDict = JSON.parsefile(pwd()*"/AngularNew.json")
    CosmoDict = JSON.parsefile(pwd()*"/Cosmology.json")
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

    Fisher = CosmoCentral.ForecastFisherαβ(ForecastContainer ,PathCentralCℓ, Path∂Cℓ,
    CosmologicalGrid, "Lensing", "Test")

    @test isapprox(39.02, CosmoCentral.EvaluateFoM(Fisher, "w0", "wa"), rtol=1e-4)
    @test isapprox(0.5670, Fisher.MarginalizedErrors["wa"], rtol=1e-4)
end