module StructuredLightFunctions

import FromFile: @from

#@from "utilities/auxiliary_functions.jl" using AuxiliaryFunctions
@from "beam_amplitude_functions.jl" using AmplitudeDistributions
#include("utilities/auxiliary_functions.jl")

#include("beam_amplitude_functions.jl")

export AmplitudeDistributions
#println(AmplitudeDistributions.Gaussian(0,0,0,1))
#println(AmplitudeDistributions.LaguerreGauss(0,0,0,1,0,0))
#println(AmplitudeDistributions.HermiteGauss(0,0,0,1,0,0))

end
