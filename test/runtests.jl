using CosmoCentral
using Test
using QuadGK
using Conda
using NumericalIntegration
Conda.add("numpy")
ENV["PYTHON"]=""
using Pkg
Pkg.build("PyCall")
using PyCall
numpy = pyimport("numpy")

params = CosmoCentral.w0waCDMParameters()
density = CosmoCentral.AnalitycalDensityStruct()
convolveddensity = CosmoCentral.ConvolvedDensityStruct()
cosmogrid  = CosmoCentral.CosmoGridStruct(zgrid=Array(LinRange(0.0, 1., 10)))

@testset "Adimensional Hubble parameter at redshift zero is equal to one" begin
    test_E_z = CosmoCentral.ComputeAdimensionalHubbleFactor(0., params)
    @test test_E_z == 1.
end

@testset "Adimensional Hubble parameter at redshift zero is equal to default value" begin
    test_H_z = CosmoCentral.ComputeHubbleFactor(0., params)
    @test test_H_z == params.H0
end

@testset "Comoving distance at redshift zero is equal to zero" begin
    test_r_z = CosmoCentral.ComputeComovingDistance(0., params)
    @test test_r_z == 0.
end

@testset "Check the normalization of density function" begin
    CosmoCentral.NormalizeAnalitycalDensityStruct(density)
    int, err = QuadGK.quadgk(x -> CosmoCentral.ComputeDensityFunction(x, density),
    density.zmin, density.zmax, rtol=1e-12)
    @test isapprox(int, density.surfacedensity, atol=1e-9)
end

@testset "Check the normalization of convolved density function" begin
    test_normalization = zeros(length(convolveddensity.zbinarray)-1)
    CosmoCentral.NormalizeConvolvedDensityStruct(convolveddensity)
    for idx in 1:length(test_normalization)
        int, err = QuadGK.quadgk(x ->
        CosmoCentral.ComputeConvolvedDensityFunction(x, idx,
        convolveddensity),
        convolveddensity.AnalitycalDensity.zmin,
        convolveddensity.AnalitycalDensity.zmax, rtol=1e-12)
        test_normalization[idx] = int
    end
    @test isapprox(test_normalization, ones(length(test_normalization)), atol=1e-12)
end

@testset "Check the computation of the convolved density function on grid" begin
    test_array = zeros(Float64,length(convolveddensity.zbinarray)-1,
    length(cosmogrid.zgrid))
    CosmoCentral.ComputeDensityFunctionConvolvedGrid(cosmogrid, convolveddensity)
    for idx_zbinarray in 1:length(convolveddensity.zbinarray)-1
        for idx_zgrid in 1:length(cosmogrid.zgrid)
            test_array[idx_zbinarray, idx_zgrid] =
            CosmoCentral.ComputeConvolvedDensityFunction(
            cosmogrid.zgrid[idx_zgrid],
            idx_zbinarray,
            convolveddensity)
        end
    end
    @test isapprox(test_array, convolveddensity.densitygridarray, atol=1e-12)
end

@testset "Check the LogSpace function against the Python equivalent" begin
    minarray = 1e-5
    maxarray = 10.
    n = 100
    becnhmark_numpy = numpy.logspace(log10(minarray), log10(maxarray), n)
    test_logspace = CosmoCentral.LogSpaced(minarray, maxarray, n)
    @test isapprox(becnhmark_numpy, test_logspace, atol=1e-12)
end
