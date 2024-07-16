using StructuredLightFunctions
using Test

@testset "StructuredLightFunctions.jl" begin
    @test AmplitudeDistributions.Gaussian(0,0,0,1) != 0
    @test AmplitudeDistributions.LaguerreGauss(0,0,0,1,0,0) != 0
    @test AmplitudeDistributions.HermiteGauss(0,0,0,1,0,0) != 0
    # Write your tests here.
end
