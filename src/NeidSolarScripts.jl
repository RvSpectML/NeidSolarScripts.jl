
#__precompile__() # this module is safe to precompile
module NeidSolarScripts
#using PyCall

#=
function __init__()
  pyimport("astropy")
end
=#

using Dates, LinearAlgebra

# For computing effects of apparent solar rotation rate
include("solar_rotation.jl")
export SolarRotation
export CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic, CalculateFWHMDifference_SolarRotation_from_loc_Equatorial
export CalculateFWHMDifference_SolarRotation_from_long_lat_alt, CalculateFWHMDifference_SolarRotation_from_obs
#export get_wiyn_loc, get_harpsn_loc

# For computing differential extinction
include("diff_solar_extinction.jl")
export DifferentialExtinction

#include("continuum_AFS.jl")
#export Continuum_AFS
include("continuum_rassine_like.jl")
export Continuum

end
