#__precompile__() # this module is safe to precompile
module NeidSolarScripts

using Dates, LinearAlgebra
#using PyCall

# For computing effects of apparent solar rotation rate
include("solar_rotation.jl")
export SolarRotation
export CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic, CalculateFWHMDifference_SolarRotation_from_loc_Equatorial
export CalculateFWHMDifference_SolarRotation_from_long_lat_alt, CalculateFWHMDifference_SolarRotation_from_obs
#export get_wiyn_loc, get_harpsn_loc

# For computing differential extinction


end
