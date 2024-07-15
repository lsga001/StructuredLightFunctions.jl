using StructuredLightFunctions
using Test

@testset "StructuredLightFunctions.jl" begin
    @test AmplitudeDistributions.LaguerreGauss(0,0,0,1,1,1)==0+0im
    @test AmplitudeDistributions.HermiteGauss(0,0,0,1,1,1)==nothing
    # Write your tests here.
end
