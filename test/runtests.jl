using StructuredLightFunctions
using Test

@testset "StructuredLightFunctions.jl" begin
    @test AmplitudeDistributions.LaguerreGauss(0,0,0,1,1,1)==0+0im
    # Write your tests here.
end
