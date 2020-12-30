using CosmoCentral
using Test

cosmo = CosmoCentral.w0waCDMParameters()

test_E_z = CosmoCentral.ComputeAdimensionalHubbleFactor(0., cosmo)
@test test_E_z == 1.

test_H_z = CosmoCentral.ComputeHubbleFactor(0., cosmo1)
@test test_H_z == 67.

test_r_z = CosmoCentral.ComputeComovingDistance(0., cosmo)
@test test_r_z == 0.
