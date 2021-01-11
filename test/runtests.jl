using CosmoCentral
#include("/home/mbonici/Desktop/CosmoCentral.jl/src/CosmoCentral.jl")
using Test
using QuadGK
using NumericalIntegration
using PyCall

numpy = pyimport("numpy")
classy = pyimport("classy")

w0waCDMCosmology = CosmoCentral.w0waCDMCosmologyStruct()
AnalitycalDensity = CosmoCentral.AnalitycalDensityStruct()
InstrumentResponse = CosmoCentral.InstrumentResponseStruct()
ConvolvedDensity = CosmoCentral.ConvolvedDensityStruct()
CosmologicalGrid  = CosmoCentral.CosmologicalGridStruct(
ZArray=Array(LinRange(0.001, 2.5, 300)))
BackgroundQuantities = CosmoCentral.BackgroundQuantitiesStruct(HZArray=
zeros(length(CosmologicalGrid.ZArray)),
rZArray=zeros(length(CosmologicalGrid.ZArray)))
PiecewiseBias = CosmoCentral.PiecewiseBiasStruct(BiasArray =
zeros(length(ConvolvedDensity.DensityNormalizationArray),
length(CosmologicalGrid.ZArray)))
GCWeightFunction = CosmoCentral.GCWeightFunctionStruct(WeightFunctionArray =
zeros(length(ConvolvedDensity.DensityNormalizationArray),
length(CosmologicalGrid.ZArray)))
AngularCoefficients = CosmoCentral.AngularCoefficientsStruct(
AngularCoefficientsArray = zeros(length(CosmologicalGrid.MultipolesArray),
length(GCWeightFunction.WeightFunctionArray[:, 1]),
length(GCWeightFunction.WeightFunctionArray[:, 1])))

@testset "Evaluation of background quantities" begin
    test_E_z = CosmoCentral.ComputeAdimensionalHubbleFactor(0., w0waCDMCosmology)
    @test test_E_z == 1.
    test_H_z = CosmoCentral.ComputeHubbleFactor(0., w0waCDMCosmology)
    @test test_H_z == w0waCDMCosmology.H0
    test_r_z = CosmoCentral.ComputeComovingDistance(0., w0waCDMCosmology)
    @test test_r_z == 0.
    test_r_array = zeros(length(CosmologicalGrid.ZArray))
    test_H_array = zeros(length(CosmologicalGrid.ZArray))
    for (idxz, zvalue) in enumerate(CosmologicalGrid.ZArray)
        test_r_array[idxz] = CosmoCentral.ComputeComovingDistance(zvalue,
        w0waCDMCosmology)
        test_H_array[idxz] = CosmoCentral.ComputeHubbleFactor(zvalue,
        w0waCDMCosmology)
    end
    CosmoCentral.ComputeBackgroundQuantitiesOverGrid(CosmologicalGrid,
    BackgroundQuantities, w0waCDMCosmology)
    @test test_r_array ==BackgroundQuantities.rZArray
    @test test_H_array ==BackgroundQuantities.HZArray
end

@testset "Check the normalization of density function" begin
    CosmoCentral.NormalizeAnalitycalDensityStruct(AnalitycalDensity)
    int, err = QuadGK.quadgk(x -> CosmoCentral.ComputeDensityFunction(x, AnalitycalDensity),
    AnalitycalDensity.ZMin, AnalitycalDensity.ZMax, rtol=1e-12)
    @test isapprox(int, AnalitycalDensity.SurfaceDensity, atol=1e-9)
end

@testset "Check the normalization of convolved density function" begin
    test_normalization = zeros(length(ConvolvedDensity.ZBinArray)-1)
    CosmoCentral.NormalizeConvolvedDensityStruct(ConvolvedDensity, AnalitycalDensity,
    InstrumentResponse, CosmologicalGrid)
    for idx in 1:length(test_normalization)
        int, err = QuadGK.quadgk(x ->
        CosmoCentral.ComputeConvolvedDensityFunction(x, idx,
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
    CosmoCentral.ComputeConvolvedDensityFunctionGrid(CosmologicalGrid, ConvolvedDensity,
    AnalitycalDensity, InstrumentResponse)
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            test_array[idx_ZBinArray, idx_ZArray] =
            CosmoCentral.ComputeConvolvedDensityFunction(
            CosmologicalGrid.ZArray[idx_ZArray],
            idx_ZBinArray,
            ConvolvedDensity,
            AnalitycalDensity,
            InstrumentResponse)
        end
    end
    @test isapprox(test_array, ConvolvedDensity.DensityGridArray, atol=1e-12)
end

@testset "Check the LogSpace function against the Python equivalent" begin
    minarray = 1e-5
    maxarray = 10.
    n = 100
    benchmark_numpy = numpy.logspace(log10(minarray), log10(maxarray), n)
    test_logspace = CosmoCentral.LogSpaced(minarray, maxarray, n)
    @test isapprox(benchmark_numpy, test_logspace, atol=1e-12)
end

@testset "Check the Binsearch function" begin
    array = [1.,2.,3.]
    @test 1 == CosmoCentral.BinSearch(1.5, array)
end

@testset "Check the Piecewise bias evaluation" begin
    z = 0.4
    test_array = zeros(length(CosmologicalGrid.ZArray))
    bias = CosmoCentral.ComputeBias(z, PiecewiseBias, ConvolvedDensity)
    @test isapprox(1.0997727037892875, bias, atol=1e-12)
    for (idxz, zvalue) in enumerate(CosmologicalGrid.ZArray)
        test_array[idxz] = CosmoCentral.ComputeBias(zvalue, PiecewiseBias,
        ConvolvedDensity)
    end
    CosmoCentral.ComputeBiasOverGrid(CosmologicalGrid, PiecewiseBias,
    ConvolvedDensity)
    @test isapprox(test_array, PiecewiseBias.BiasArray[1,:], atol=1e-12)
end

@testset "Check the Weight function evaluation" begin
    test_array = zeros(length(ConvolvedDensity.ZBinArray)-1,
    length(CosmologicalGrid.ZArray))
    for (idxz, zvalue) in enumerate(CosmologicalGrid.ZArray)
        for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
            test_array[idx_ZBinArray, idxz] =
            CosmoCentral.ComputeWeightFunction(zvalue, idx_ZBinArray,
        ConvolvedDensity, AnalitycalDensity,
        InstrumentResponse, w0waCDMCosmology,
        PiecewiseBias, GCWeightFunction)
        end
    end
CosmoCentral.ComputeWeightFunctionOverGrid(GCWeightFunction, AnalitycalDensity, InstrumentResponse, ConvolvedDensity, PiecewiseBias, CosmologicalGrid, BackgroundQuantities, w0waCDMCosmology)
@test isapprox(test_array, GCWeightFunction.WeightFunctionArray, atol=1e-9)
end


@testset "Check the Power Spectrum evaluated over the Limber Grid" begin
    classyParams = CosmoCentral.Initializeclassy(w0waCDMCosmology)
    PowerSpectrum = CosmoCentral.PowerSpectrumStruct(PowerSpectrumLinArray =
    zeros(length(CosmologicalGrid.KArray), length(CosmologicalGrid.ZArray)),
    PowerSpectrumNonlinArray = zeros(length(CosmologicalGrid.KArray),
    length(CosmologicalGrid.ZArray)),
    InterpolatedPowerSpectrum = zeros(length(CosmologicalGrid.MultipolesArray),
    length(CosmologicalGrid.ZArray)))
    CosmoCentral.ComputeLimberArray(CosmologicalGrid, BackgroundQuantities)
    test_k_limber = (CosmologicalGrid.MultipolesArray[1]+0.5) /
    BackgroundQuantities.rZArray[1]
    @test test_k_limber == CosmologicalGrid.KLimberArray[1, 1]
    test_Omega_cdm = w0waCDMCosmology.ΩM-w0waCDMCosmology.ΩB-
    w0waCDMCosmology.Mν/(93.14*(w0waCDMCosmology.H0/100)^2)
    @test test_Omega_cdm == classyParams.classyParamsDict["Omega_cdm"]
    #CosmoCentral.EvaluatePowerSpectrum(classyParams, CosmologicalGrid,
    #PowerSpectrum)
    #CosmoCentral.InterpolateAndEvaluatePowerSpectrum(CosmologicalGrid,
    #BackgroundQuantities, PowerSpectrum)
    #cosmo = classy.Class()
    #cosmo.set(classyParams.classyParamsDict)
    #cosmo.compute()
    #println(cosmo.pk(test_k_limber, CosmologicalGrid.ZArray[1]))
    #test_power_spectrum = cosmo.pk(test_k_limber, CosmologicalGrid.ZArray[1])
    #println(test_power_spectrum)
    #@test isapprox(test_power_spectrum,
    #PowerSpectrum.InterpolatedPowerSpectrum[1, 1], atol=1e-7)
end
