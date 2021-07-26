using CosmoCentral
using Test
using QuadGK
using NumericalIntegration
using PyCall
using Interpolations

numpy = pyimport("numpy")

w0waCDMCosmology = CosmoCentral.Flatw0waCDMCosmology()
AnalitycalDensity = CosmoCentral.AnalitycalDensity()
InstrumentResponse = CosmoCentral.InstrumentResponse()
ConvolvedDensity = CosmoCentral.ConvolvedDensity()
CosmologicalGrid  = CosmoCentral.CosmologicalGrid(
ZArray=Array(LinRange(0.001, 4.0, 500)))
BackgroundQuantities = CosmoCentral.BackgroundQuantities(HZArray=
zeros(length(CosmologicalGrid.ZArray)),
rZArray=zeros(length(CosmologicalGrid.ZArray)))
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
input_path_pmm = pwd()*"/p_mm"
input_path_Cℓ = pwd()*"/cl"
CosmoCentral.classy.Class()

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
    @test test_r_array ==BackgroundQuantities.rZArray
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
    BackgroundQuantitiesLoaded.rZArray[1]
    @test test_k_limber == CosmologicalGrid.KLimberArray[1, 1]
    test_Omega_cdm = w0waCDMCosmology.ΩM-w0waCDMCosmology.ΩB-
    w0waCDMCosmology.Mν/(93.14*(w0waCDMCosmology.H0/100)^2)
    @test test_Omega_cdm == classyParams.classyParamsDict["Omega_cdm"]
    x = LinRange(log10(first(CosmologicalGrid.KArray)),
    log10(last(CosmologicalGrid.KArray)), length(CosmologicalGrid.KArray))
    y = LinRange(first(CosmologicalGrid.ZArray), last(CosmologicalGrid.ZArray),
    length(CosmologicalGrid.ZArray))
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
    @test isapprox(BackgroundQuantitiesLoaded.rZArray,
    NewBackgroundQuantitiesLoaded.rZArray, rtol=1e-6)
    PowerSpectrumDierckx, BackgroundQuantitiesLoaded, CosmologicalGrid =
    CosmoCentral.ReadPowerSpectrumBackground(input_path_pmm, ℓBinCenters, Array(ones(2990)))
    CosmoCentral.InterpolatePowerSpectrumLimberGrid!(CosmologicalGrid,
    BackgroundQuantitiesLoaded, PowerSpectrumDierckx,
    CosmoCentral.RectBivSplineDierckx())
    @test isapprox(PowerSpectrumDierckx.InterpolatedPowerSpectrum[1, 1],
    PowerSpectrum.InterpolatedPowerSpectrum[1, 1], rtol=1e-2)
end

@testset "Test Angular coefficients evaluation" begin
    CℓLoaded = CosmoCentral.ReadCℓ(input_path_Cℓ, "PhotometricGalaxy_PhotometricGalaxy")
    ℓBinCenters = Array(LinRange(10.5, 2999.5, 2990))
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
