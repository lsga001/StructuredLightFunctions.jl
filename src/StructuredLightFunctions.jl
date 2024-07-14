module StructuredLightFunctions

export LaguerreGauss

using BeamAmplitudeFunctions

println(AmplitudeDistributions.LaguerreGauss(0,0,0,1,1,1))
println(AmplitudeDistributions.HermiteGauss(0,0,0,1,1,1))

end
