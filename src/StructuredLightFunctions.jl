module StructuredLightFunctions

import FromFile: @from

#@from "utilities/auxiliary_functions.jl" using AuxiliaryFunctions
@from "beam_amplitude_functions.jl" using AmplitudeDistributions
#include("utilities/auxiliary_functions.jl")

#include("beam_amplitude_functions.jl")

export AmplitudeDistributions

end
