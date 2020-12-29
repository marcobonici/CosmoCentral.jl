import CosmoCentral
using Test

ciccio = CosmoCentral.ComputeHubbleFactor(0.)

@test ciccio == 0.
@test ciccio == 0.
