using CosmoCentral
using Test

test_E_z = CosmoCentral.ComputeAdimensionalHubbleFactor(0.)
@test test_E_z == 1.

test_H_z = CosmoCentral.ComputeHubbleFactor(0.)
@test test_H_z == 67.

test_r_z = CosmoCentral.ComputeComovingDistance(0.)
@test test_r_z == 0.
