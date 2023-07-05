
#__precompile__() # this module is safe to precompile
module NeidSolarScripts
#using PyCall

#=
function __init__()
  pyimport("astropy")
end
=#

using Dates, LinearAlgebra

# For reading data from pyrheliometer files
include("pyroheliometer.jl")
export Pyroheliometer

#=
# For computing effects of apparent solar rotation rate
include("solar_rotation.jl")
export SolarRotation
export CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic, CalculateFWHMDifference_SolarRotation_from_loc_Equatorial
export CalculateFWHMDifference_SolarRotation_from_long_lat_alt, CalculateFWHMDifference_SolarRotation_from_obs
#export get_wiyn_loc, get_harpsn_loc

# For computing differential extinction
include("diff_solar_extinction.jl")
export DifferentialExtinction
=#

#include("continuum_AFS.jl")
#export Continuum_AFS
include("continuum_rassine_like.jl")
export Continuum


include("arg_parse.jl")
export parse_commandline_make_manifest_solar, parse_commandline_make_pyrheliometer_daily, parse_commandline_verify_downloads
export parse_commandline_calc_order_ccfs 
export parse_commandline_daily_report, parse_commandline_combine_daily_reports, parse_commandline_combine_daily_rvs

end
