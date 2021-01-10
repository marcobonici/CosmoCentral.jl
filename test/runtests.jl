using CosmoCentral
using Test
using QuadGK
using NumericalIntegration
using PyCall
numpy = pyimport("numpy")

w0waCDMCosmology = CosmoCentral.w0waCDMCosmologyStruct()
AnalitycalDensity = CosmoCentral.AnalitycalDensityStruct()
instrumentresponse = CosmoCentral.InstrumentResponseStruct()
ConvolvedDensity = CosmoCentral.ConvolvedDensityStruct()
CosmologicalGrid  = CosmoCentral.CosmologicalGridStruct(ZArray=Array(LinRange(0.2, 1., 100)))
PiecewiseBias = CosmoCentral.PiecewiseBiasStruct(BiasArray = zeros(length(ConvolvedDensity.DensityNormalizationArray), length(CosmologicalGrid.ZArray)))

@testset "Adimensional Hubble parameter at redshift zero is equal to one" begin
    test_E_z = CosmoCentral.ComputeAdimensionalHubbleFactor(0., w0waCDMCosmology)
    @test test_E_z == 1.
end

@testset "Adimensional Hubble parameter at redshift zero is equal to default value" begin
    test_H_z = CosmoCentral.ComputeHubbleFactor(0., w0waCDMCosmology)
    @test test_H_z == w0waCDMCosmology.H0
end

@testset "Comoving distance at redshift zero is equal to zero" begin
    test_r_z = CosmoCentral.ComputeComovingDistance(0., w0waCDMCosmology)
    @test test_r_z == 0.
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
    instrumentresponse)
    for idx in 1:length(test_normalization)
        int, err = QuadGK.quadgk(x ->
        CosmoCentral.ComputeConvolvedDensityFunction(x, idx,
        ConvolvedDensity, AnalitycalDensity, instrumentresponse),
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
    AnalitycalDensity, instrumentresponse)
    for idx_ZBinArray in 1:length(ConvolvedDensity.ZBinArray)-1
        for idx_ZArray in 1:length(CosmologicalGrid.ZArray)
            test_array[idx_ZBinArray, idx_ZArray] =
            CosmoCentral.ComputeConvolvedDensityFunction(
            CosmologicalGrid.ZArray[idx_ZArray],
            idx_ZBinArray,
            ConvolvedDensity,
            AnalitycalDensity,
            instrumentresponse)
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
