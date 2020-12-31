using CosmoCentral
using Test

params = CosmoCentral.w0waCDMParameters()
density = CosmoCentral.AnalitycalDensityStruct()

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
    test_density = CosmoCentral.NormalizeAnalitycalDensityStruct(density)
    normalization, err = QuadGK.quadgk(x -> 1 /
    ComputeDensityFunction(x, densityparameters), densityparameters.zmin,
    densityparameters.zmax, rtol=1e-12)*densityparameters.surfacedensity
    @test normalization == test_density.normalization
end
