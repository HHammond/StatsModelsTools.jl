module StatsModelsTools

using StatsBase
using StatsModels

include("rms_splines.jl")
include("polynomials.jl")

using .RMSSplines
using .Polynomials

export rcs, rms_generate_knots, rms_quantiles_for_knots
export poly

end # module
