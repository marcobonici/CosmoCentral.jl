import CosmoCentral
using Test

ciccio = CosmoCentral.ComputeAdimensionalHubbleFactor(0.)

@test ciccio == 1.

ciccio = CosmoCentral.ComputeHubbleFactor(0.)
@test ciccio == 67.
println("Everything work!")
