verbose = false
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectMLBase
 using EchelleInstruments#, EchelleInstruments.NEID
 using RvSpectML
 if verbose println("# Loading NeidSolarScripts")    end
 using SunAsAStar
 using NeidSolarScripts
 #using ArgParse

println("# Parsing arguments...")
   lsf_width = 3.0e3
#=
function parse_commandline()
     s = ArgParseSettings( description = "Calculate order CCFs from NEID L2 FITS files.")
     #import_settings!(s, s_files_only, args_only=false)
     add_arg_group!(s, "Files", :argg_files)
     @add_arg_table! s begin
         "manifest"
             help = "Manifest file (CVS) containing FITS files to analyze.\nExpects columns named Filename, bjd, target, and used for continuum continuum_filename."
             arg_type = String
             default = "manifest.csv"
             #required = true
         "output"
             help = "Filename for output CCFs (jld2)"
             arg_type = String
             default = "daily_ccfs.jld2"
         "--param"
             help = "Parameters file"
             arg_type = String
         "--overwrite"
            help = "Specify it's ok to overwrite the output file."
            action = :store_true
         "--verbose"
            help = "Verbose logging."
            action = :store_true
      end
      add_arg_group!(s, "CCF parameters", :argg_ccf_param)
      @add_arg_table! s begin
        "--orders_to_use"
           help = "First and last order _index_ to compute CCFs for"
           nargs = 2
           arg_type = Int64
           #default = [ min_order(NEID2D()), max_order(NEID2D()) ]
           #default = [ first(orders_to_use_default(NEID2D())), last(orders_to_use_default(NEID2D())) ]
        "--mask_scale_factor"
            help = "Specify CCF mask width scale as multiple of NEID default v width " * string(default_ccf_mask_v_width(NEID2D()))
            arg_type = Float64
            #default = round(lsf_width/default_ccf_mask_v_width(NEID2D()), sigdigits=3)
        "--line_width_50_default"
            help = "Specify default line full width half maximum"
            arg_type = Float64
            #default = 7.9e3 # m/2
        "--variable_mask_scale"
            help = "Vary width of mask to compensate for solar rotation."
            action = :store_true
        "--ccf_mid_vel"
            help = "Middle of velocity range to use for CCFs."
            arg_type = Float64
            #default = 0.0
        "--range_no_mask_change"
            help = "Avoid lines that would result in a mask change due to BC using line width times this factor (TODO: Verify definition)"
            arg_type = Float64
            #default = 5.0
        "--v_step"
            help = "Specify v step for computing CCF "
            arg_type = Float64
            #default = 155.0 # m/2
      end
      add_arg_group!(s, "Line list parameters", :argg_line_list_param)
      @add_arg_table! s begin
         "--line_list_filename"
             help = "Line list filename (input)"
             arg_type = String
             #default = joinpath(pkgdir(NeidSolarScripts),"data","solar_line_list_espresso.csv")
         "--line_list_output_filename"
             help = "Line list filename (output)"
             arg_type = String
         "--recompute_line_weights"
             help = "Force recalculation of SNR-based line weight factors."
             action = :store_true
      end
      add_arg_group!(s, "Continuum normalization parameters", :argg_continuum_param)
      @add_arg_table! s begin
         "--sed_filename"
             help = "Filename for input SED to normalize by"
             arg_type = String
             #default = "/home/eford/Code/RvSpectMLEcoSystem/NeidSolarScripts/data/neidMaster_HR_SmoothLampSED_20210101.fits"
         "--continuum_poly_half_width"
             help = "Half width for window to use for smoothing prior to polynomial fit to continuum."
             arg_type = Int64
             default = 50
         "--quantile_fit_continuum"
             help = "Quantile of rolling window to use prior to polynomial fit to continuum."
             arg_type = Float64
             default = 0.9
         "--order_poly_continuum"
             help = "Order of polynomial to fit after dividing by SED proivded."
             arg_type = Int64
             default = 4
        "--anchors_filename"
             help = "Filename for anchor locations to use in computing continuum to normalize by."
             arg_type = String
             #default = "/home/eford/Code/RvSpectMLEcoSystem/NeidSolarScripts/data/neidMaster_HR_SmoothLampSED_20210101.fits"
        "--anchors_filename_output"
             help = "Filename to write anchor locations to use in computing continuum to normalize by."
             arg_type = String
         "--smoothing_half_width"
             help = "Half width for window to use for smoothing prior to findind local maximum."
             arg_type = Int64
             default = 6
         "--stretch_factor"
             help = "Stretch factor to use for scaling flux in distance calculation for rolling pin continuum normalization"
             arg_type = Float64
             default = 5.0
         "--merging_threshold"
             help = "Threshold (in λ units) for merging nearby continuum anchors."
             arg_type = Float64
             default = 0.25
         "--fwhm_continuum"
            help = "Full width half-max to use for performing continuum normalization."
            arg_type = Float64
            default = Continuum.fwhm_sol/1000
         "--min_rollingpin_r"
            help = "Minimum rolling pin radius for continuum normalization in uints of FWHM."
            arg_type = Float64
            default = 100.0
         "--nu_continuum"
            help = "Exponent for rolling pin radius in continuum normalization"
            arg_type = Float64
            default = 1.0
         "--apply_continuum_normalization"
            help = "Apply continuum normalization."
            action = :store_true
        "--continuum_normalization_individually"
            help = "Calculate continuum normalization for each file rather than averaged spectra."
            action = :store_true
        "--save_continuum_normalization"
            help = "NOT IMPLEMENTED.  Save continuum normalization to file."
            action = :store_true
        "--recompute_continuum_normalization"
            help = "CURRENTLY ONLY OPTION.  Recompute continuum normalization even if continuum file exists."
            action = :store_true
      end
      add_arg_group!(s, "Filter manifest for ", :argg_filter_param)
      @add_arg_table! s begin
         "--target"
            help = "Target field"
            arg_type = String
            default = "Sun"
         "--datestr"
            help = "Filenames that contain date string."
            arg_type = String
        "--max_solar_hour_angle"
           help = "Maximum absolute value of solar hour angle."
           arg_type = Float64
           default = 3.12
        "--max_solar_hour_angle_clean"
           help = "Maximum absolute value of solar hour angle for clean spectrum."
           arg_type = Float64
           default = 2.0
         "--max_airmass"
            help = "Maximum airmass."
            arg_type = Float64
         "--max_airmass_clean"
            help = "Maximum airmass for clean spectrum."
            arg_type = Float64
            default = 2.0
         "--min_expmeter"
            help = "Minimum mean exposure meter flux."
            arg_type = Float64
            default = 6e4
         "--min_expmeter_clean"
            help = "Minimum mean exposure meter flux for clean spectrum."
            arg_type = Float64
            default = 1e5
         "--min_pyrhelio"
            help = "Minimum mean pyrheliometer flux."
            arg_type = Float64
            default = 10^2.95
         "--min_pyrhelio_clean"
            help = "Minimum mean pyrheliometer flux for clean spectrum."
            arg_type = Float64
            default = 10^2.95
         "--max_expmeter_rms_frac"
            help = "Maximum fractional RMS exposure meter flux."
            arg_type = Float64
            default = 0.003
         "--max_expmeter_rms_frac_clean"
            help = "Maximum fractional RMS exposure meter flux for clean spectrum."
            arg_type = Float64
            default = 0.003
         "--max_pyrhelio_rms_frac"
            help = "Maximum fractional RMS pyrheliometer flux."
            arg_type = Float64
            default = 0.0035
         "--max_pyrhelio_rms_frac_clean"
            help = "Maximum fractional RMS pyrheliometer flux for clean spectrum."
            arg_type = Float64
            default = 0.0035
         "--min_expmeter_to_pyrhelio"
            help = "Minimum mean exposure meter flux relative to pyrheliometer flux."
            arg_type = Float64
            default = 0.0
        "--min_expmeter_to_pyrhelio_clean"
            help = "Minimum mean exposure meter flux relative to pyrheliometer flux for clean spectrum."
            arg_type = Float64
            default = 0.0 # 150 from DRP v1.1
#=
         "--min_snr_factor"
            help = "Minimum SNR relative to max_snr."
            arg_type = Float64
         "--min_snr_factor_clean"
            help = "Minimum SNR relative to max_snr for clean spectrum."
            arg_type = Float64
            default = 0.5
         "--max_snr"
            help = "Specify max_snr manually."
            arg_type = Float64
=#
         "--start_time"
            help = "Specify daily start time for CCFs to be processed (HH MM)"
            nargs = 2
            arg_type = Int64
            default = [0, 0]
            #default = [18, 30]
        "--stop_time"
           help = "Specify daily stop time for CCFs to be processed (HH MM)"
           nargs = 2
           arg_type = Int64
           #default = [ min_order(NEID2D()), max_order(NEID2D()) ]
           #default = [22, 30]
           default = [23, 59]
         "--start_time_clean"
            help = "Specify daily start time for CCFs for clean spectrum"
            nargs = 2
            arg_type = Int64
            default = [17, 30]
        "--stop_time_clean"
           help = "Specify daily stop time for CCFs for clean spectrum"
           nargs = 2
           arg_type = Int64
           default = [22, 12]
         "--max_num_spectra"
            help = "Maximum number of spectra to process."
            arg_type = Int
            default = 300  # Enough for one day of NEID
         "--max_num_spectra_clean"
            help = "Maximum number of spectra to include in clean spectrum."
            arg_type = Int
            default = 130  # based on 120 minutes of integration time from 55s exposures
            #default = 65  # based on 60 minutes of integration time from 55s exposures
      end

     return parse_args(s)
 end
 args = parse_commandline()
=#
 args = parse_commandline_calc_order_ccfs()

verbose = haskey(args,"verbose") ? args["verbose"] : verbose

if verbose println(now()) end
 #println("# Loading other packages 1/2")

if verbose   println("# Loading other packages 2/2")    end
 using CSV, DataFrames, Query, Dates
 using JLD2, FileIO, MD5 #, SHA
 using StatsBase, Statistics, NaNMath

function extract_time_from_filename(fn::AbstractString)
    m = match(r"neidL\d_(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})(\d{2})\.fits$",basename(fn))
    if length(m)<6 
        @warn("Can't extract time from " * fn)
        return nothing
    end
    return Time(parse(Int,m[4]),parse(Int,m[5]),parse(Int,m[6]))
end
    
 # Filename arguments
 @assert isfile(args["manifest"]) || islink(args["manifest"])
 manifest_filename = args["manifest"]
 #args["overwrite"] = true   # TODO: Comment after done testing
 if isfile(args["output"]) && !args["overwrite"]
    error("# Output file " * args["output"] * " already exists (size = " * string(filesize(args["output"])) * " ).")
 end
 @assert !isfile(args["output"]) || args["overwrite"] == true
 @assert match(r"\.jld2$",args["output"]) != nothing
 daily_ccf_filename = args["output"]
 touch(daily_ccf_filename)

 file_hashes = Dict{String,String}()
 start_processing_time = now()

 # Include parameters file if specified
 if args["param"] != nothing
     @warn "Haven't implemented parameters file for this script yet."
    if isfile(args["param"]) || islink(args["param"])
        println("# Including parameter code from ", args["param"], ".")
        include(args["param"])
        file_hashes["param.in"] = bytes2hex(args["param"])
    end
 end
 # Overide values from param file if on commandline.  If missing, set defaults
 if args["line_list_filename"] != nothing
    line_list_filename = args["line_list_filename"]
 elseif !@isdefined line_list_filename
     #line_list_filename = joinpath(pkgdir(NeidSolarScripts),"data","solar_line_list_espresso.csv")
     line_list_filename = joinpath(pkgdir(NeidSolarScripts),"data","espresso+neid_mask_97_to_108.mas")
 end
 @assert isfile(line_list_filename) || islink(line_list_filename)
 if args["sed_filename"] != nothing
    @warn("DRP v1.1 now provides blaze in L2 file.  This script has not been updated to use explicit an SED model.") 
    sed_filename = args["sed_filename"]
 #elseif !@isdefined sed_filename
 #     sed_filename = joinpath("/home/eford/Code/RvSpectMLEcoSystem/NeidSolarScripts/data","neidMaster_HR_SmoothLampSED_20210101.fits")
 end
 @assert !@isdefined(sed_filename) || isfile(sed_filename) || islink(sed_filename)
 @assert !@isdefined(anchors_filename) || isfile(anchors_filename) || islink(anchors_filename)

 # Extract arguments with non-standard types
 if args["orders_to_use"] != nothing && length(args["orders_to_use"]) == 2
     orders_to_use = first(args["orders_to_use"]):last(args["orders_to_use"])
 elseif !@isdefined orders_to_use
     orders_to_use = first(orders_to_use_default(NEID2D())):last(orders_to_use_default(NEID2D()))
 end
 @assert min_order(NEID2D()) <= first(orders_to_use) <= max_order(NEID2D())
 @assert min_order(NEID2D()) <= last(orders_to_use) <= max_order(NEID2D())


 if args["ccf_mid_vel"] != nothing
     ccf_mid_velocity = args["ccf_mid_vel"]
 elseif !@isdefined ccf_mid_velocity
     ccf_mid_velocity = 0.0
 end
 @assert -5e5 < ccf_mid_velocity < 5e5  # arbitrary ball park restrictions to catch obvious typeos

 if args["v_step"] != nothing
     v_step = args["v_step"]
 elseif !@isdefined v_step
     v_step = 150.0
 end
 @assert 50 <= v_step <= 2000 # m/s

 if args["mask_scale_factor"] != nothing
     mask_scale_factor = args["mask_scale_factor"]
 elseif !@isdefined mask_scale_factor
     mask_scale_factor = round(lsf_width/default_ccf_mask_v_width(NEID2D()), sigdigits=3)
 end
 @assert 1 <= mask_scale_factor <= 64 # arbitrary big number to catch obvious typeos

 if args["line_width_50_default"] != nothing
     line_width_50_default = args["line_width_50_default"]
 elseif !@isdefined line_width_50_default
     line_width_50_default = 7.9e3 # m/s
 end
 @assert 2e3 <= line_width_50_default <= 50e3   # arbitrary big number to catch obvious typeos

 if args["range_no_mask_change"] != nothing
     range_no_mask_change = args["range_no_mask_change"]
 elseif !@isdefined range_no_mask_change
     range_no_mask_change = 5.0 # m/s
 end
 @assert 0 <= range_no_mask_change <= 10
 start_time = args["start_time"] != nothing && length(args["start_time"]) == 2 ? Time(args["start_time"][1], args["start_time"][2]) : Time(0,0)
 stop_time =  args["stop_time"] != nothing && length(args["stop_time"]) == 2 ? Time(args["stop_time"][1],  args["stop_time"][2]) : Time(12,59,59)
 start_time_clean = args["start_time_clean"] != nothing && length(args["start_time_clean"]) == 2 ? Time(args["start_time_clean"][1], args["start_time_clean"][2]) : Time(0,0)
 stop_time_clean =  args["stop_time_clean"] != nothing && length(args["stop_time_clean"]) == 2 ? Time(args["stop_time_clean"][1],  args["stop_time_clean"][2]) : Time(12,59,59)

 @assert 3 <= args["continuum_poly_half_width"] <= 200 # arbitrary limits for now
 @assert 0.5 <= args["quantile_fit_continuum"] <= 1-1/(2*args["smoothing_half_width"])
 @assert 0 <= args["order_poly_continuum"] <= 6 # arbitrary limits for now
 @assert 3 <= args["smoothing_half_width"] <= 200 # arbitrary limits for now
 @assert 1.0 <= args["stretch_factor"] <= 100    # arbitrary limits for now
 @assert 0.1 <= args["merging_threshold"] <= 2.0 # arbitrary limits for now
 @assert 1.0 <= args["fwhm_continuum"] <= 30.0   # km/s  arbitrary limits for now
 @assert 10 <= args["min_rollingpin_r"] <= 1000 # arbitrary limits for now
 @assert 0.7 <= args["nu_continuum"] <= 1.3  # Recommendations from Cretignier et al.
 #args["apply_continuum_normalization"] = true   # TODO Remove after done testing
 #args["continuum_normalization_individually"] = false #  TODO Remove after done testing
 #args["recompute_line_weights"] = true  # TODO Remove after done testing

if filesize(manifest_filename) == 0
   println("# Manifest filesize = 0 implies error generating manifest: ", manifest_filename, ".")
   exit(0)
end

if verbose println("# Reading manifest of files to process.")  end
  df_files  = CSV.read(manifest_filename, DataFrame)
    @assert size(df_files,1) >= 1
    @assert hasproperty(df_files,:Filename)
    @assert hasproperty(df_files,:target)
    @assert hasproperty(df_files,:bjd)
    @assert hasproperty(df_files,:ssbz)
    @assert hasproperty(df_files,:exptime)
    @assert hasproperty(df_files,:airmass)
    if args["target"] == "Sun" || args["target"] == "Solar"
      @assert hasproperty(df_files,:alt_sun)  # TODO: Compute if not avaliable?
      @assert hasproperty(df_files,:Δfwhm²)   # TODO: Compute if not avaliable?
      if !hasproperty(df_files,:hour_angle)
         df_files[!,"hour_angle"] = SolarRotation.get_solar_hour_angle(df_files.bjd,obs=:WIYN)
      end
      #if eltype(df_files[!,:order_snrs]) == String   # TODO: Compute if not avaliable?
    end
    #=
    @assert hasproperty(df_files,:order_snrs)
        df_files[!,:order_snrs] = map(i->parse.(Float64,split(df_files[i,:order_snrs][2:end-1],',')),1:size(df_files,1))
    end
    @assert eltype(df_files[!,:order_snrs]) == Vector{Float64}
    =#

  #=
  if  args["max_snr"] != nothing
      max_snr = args["max_snr"]
  elseif !@isdefined max_snr
      df_files_tmp = df_files |>
        @filter( args["target"] == nothing || _.target == args["target"] ) |> DataFrame
      max_snr = maximum(NaNMath.sum.(df_files_tmp.order_snrs))
  end
  @assert 0 < max_snr < Inf
  =#

  @assert args["max_airmass"] == nothing || 1 < args["max_airmass"] <= 10  # not tested beyond 10
  #@assert args["min_snr_factor"] == nothing || 0 <= args["min_snr_factor"] < 1
  @assert args["max_solar_hour_angle"] == nothing || 0 < args["max_solar_hour_angle"] <= 6
 
 min_drp_minor_version = VersionNumber(1,2,0)
 max_drp_minor_version = VersionNumber(1,2,0)


  df_files_use = df_files |>
    @filter( args["target"] == nothing || _.target == args["target"] ) |>
    DataFrame
  df_files_use = df_files_use |>   
    @filter( args["datestr"] == nothing || occursin(args["datestr"],_.Filename) ) |>
    DataFrame
  df_files_use = df_files_use |>   
    @filter( _.drpversion != "" && (min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version )) |> 
    DataFrame
  df_files_use = df_files_use |>   
    @filter( args["max_solar_hour_angle"] == nothing || abs(_.hour_angle) <= args["max_solar_hour_angle"] ) |>
    DataFrame
  df_files_use = df_files_use |>   
    @filter( args["max_airmass"] == nothing || _.airmass <= args["max_airmass"] ) |>
    DataFrame
  df_files_use = df_files_use |>   
    @filter( args["min_expmeter"] == nothing || _.expmeter_mean >= args["min_expmeter"] ) |>
    DataFrame
  if !hasproperty(df_files_use,:mean_pyroflux) || !hasproperty(df_files_use,:rms_pyroflux) || any(ismissing.(df_files_use.mean_pyroflux)) || any(ismissing.(df_files_use.rms_pyroflux))
     @error "Manifest file doesn't include valid mean_pyroflux and/or rms_pyroflux." manifest_filename 
  end
  df_files_use = df_files_use |>   
    @filter( args["min_pyrhelio"] == nothing || _.mean_pyroflux >= args["min_pyrhelio"] ) |>
    DataFrame
  df_files_use = df_files_use |>   
    @filter( args["max_expmeter_rms_frac"] == nothing || _.expmeter_rms <= args["max_expmeter_rms_frac"]*_.expmeter_mean  ) |>
    DataFrame
  df_files_use = df_files_use |>   
    @filter( args["max_pyrhelio_rms_frac"] == nothing || _.rms_pyroflux <= args["max_pyrhelio_rms_frac"]*_.mean_pyroflux  ) |>
    DataFrame
 df_files_use = df_files_use |>   
    @filter( args["min_expmeter_to_pyrhelio"] == nothing || _.expmeter_mean >= args["min_expmeter_to_pyrhelio"]*_.mean_pyroflux  ) |>
    DataFrame
 df_files_use = df_files_use |>   
    @filter( args["start_time"] == nothing || extract_time_from_filename(_.Filename) >= start_time ) |>
    DataFrame
  df_files_use = df_files_use |>   
    @filter( args["stop_time"] == nothing || extract_time_from_filename(_.Filename) <= stop_time ) |> # TODO for other instruments may need to deal wtih cross end of 24 UTC
    DataFrame
  df_files_use = df_files_use |>   
    @filter( args["max_airmass"] == nothing || _.airmass <= args["max_airmass"] ) |>
    DataFrame
  if hasproperty(df_files_use,:dq1level)
    df_files_use = df_files_use |>   
       @filter( mod(_.dq1level,4) <2  ) |>
       DataFrame
  end
  
  #=
 df_files_use = df_files_use |>   
    @filter( _.drpversion != "" && (min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version )) |> 
    DataFrame
 if !hasproperty(df_files_use,:mean_pyroflux) || !hasproperty(df_files_use,:rms_pyroflux) || any(ismissing.(df_files_use.mean_pyroflux)) || any(ismissing.(df_files_use.rms_pyroflux))
     @error "Manifest file doesn't include valid mean_pyroflux and/or rms_pyroflux." manifest_filename 
  end  
 df_files_use = df_files_use |>   
    @filter( _.rms_pyroflux <= 0.0035* _.mean_pyroflux ) |> 
    DataFrame
  df_files_use = df_files_use |>   
    @filter( _.expmeter_mean >= 6e4 ) |> 
    DataFrame
 =#
  df_files_use = df_files_use |>   
    @filter( _.driftfun == "dailymodel0" ) |>
    DataFrame
  df_files_use = df_files_use |>   
    @orderby(_.bjd) |>
    @take(args["max_num_spectra"] ) |>
   DataFrame
  println("# Found ", size(df_files_use,1), " files of ",  size(df_files,1), " to process.")
  if !(size(df_files_use,1) >= 1)
    println("# Found no files passing requirements.\n")
    exit(0)
  end
  @assert size(df_files_use,1) >= 1

#max_drp_minor_version = Base.thisminor(maximum(VersionNumber.(df_files_use.drpversion)))

#=
df_files_cleanest = df_files_use |>
    @filter( min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version ) |> 
    @filter( _.airmass <= args["max_airmass_clean"] ) |>
    @filter( abs(_.hour_angle) <= args["max_solar_hour_angle_clean"] ) |>
    DataFrame
  println("# Found ", size(df_files_cleanest,1), " files considered clean 1.")

df_files_cleanest = df_files_use |>
    @filter( min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version ) |> 
    @filter( _.airmass <= args["max_airmass_clean"] ) |>
    @filter( abs(_.hour_angle) <= args["max_solar_hour_angle_clean"] ) |>
    @filter( _.expmeter_mean >= args["min_expmeter_clean"] ) |> 
    @filter( _.mean_pyroflux >= args["min_pyrhelio_clean"] ) |> 
    @take(args["max_num_spectra_clean"] ) |>
    DataFrame
  println("# Found ", size(df_files_cleanest,1), " files considered clean 2.")

df_files_cleanest = df_files_use |>
    @filter( min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version ) |> 
    @filter( _.airmass <= args["max_airmass_clean"] ) |>
    @filter( abs(_.hour_angle) <= args["max_solar_hour_angle_clean"] ) |>
    @filter( _.expmeter_mean >= args["min_expmeter_clean"] ) |> 
    @filter( _.mean_pyroflux >= args["min_pyrhelio_clean"] ) |> 
    @filter( _.expmeter_rms <= args["max_expmeter_rms_frac_clean"]*_.expmeter_mean ) |> 
    @filter( _.rms_pyroflux <= args["max_pyrhelio_rms_frac_clean"]*_.mean_pyroflux  ) |> 
    DataFrame
  println("# Found ", size(df_files_cleanest,1), " files considered clean. 3")
=#

println("# Extracting time (", extract_time_from_filename(df_files_use.Filename[1]), ") from first filename (", df_files_use.Filename[1], ")")
println("# start_time_clean = ", start_time_clean)
println("# stop_time_clean = ", stop_time_clean)
for fn in df_files_use.Filename
    println("# ", start_time_clean<=extract_time_from_filename(fn)<=stop_time_clean, " time (", extract_time_from_filename(fn), ") from filename (", fn, ")")
end

df_files_cleanest = df_files_use |>
    @filter( min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version ) |> 
    @filter( _.airmass <= args["max_airmass_clean"] ) |>
    @filter( abs(_.hour_angle) <= args["max_solar_hour_angle_clean"] ) |>
    DataFrame
  println("# Found ", size(df_files_cleanest,1), " files considered clean 0.")

df_files_cleanest = df_files_use |>
    @filter( min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version ) |> 
    @filter( _.airmass <= args["max_airmass_clean"] ) |>
    @filter( abs(_.hour_angle) <= args["max_solar_hour_angle_clean"] ) |>
    @filter( _.expmeter_mean >= args["min_expmeter_clean"] ) |> 
    @filter( _.mean_pyroflux >= args["min_pyrhelio_clean"] ) |> 
    DataFrame
  println("# Found ", size(df_files_cleanest,1), " files considered clean 1.")

df_files_cleanest = df_files_use |>
    @filter( min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version ) |> 
    @filter( _.airmass <= args["max_airmass_clean"] ) |>
    @filter( abs(_.hour_angle) <= args["max_solar_hour_angle_clean"] ) |>
    @filter( _.expmeter_mean >= args["min_expmeter_clean"] ) |> 
    @filter( _.mean_pyroflux >= args["min_pyrhelio_clean"] ) |> 
    @filter( _.expmeter_rms <= args["max_expmeter_rms_frac_clean"]*_.expmeter_mean ) |> 
    @filter( _.rms_pyroflux <= args["max_pyrhelio_rms_frac_clean"]*_.mean_pyroflux  ) |> 
    DataFrame
  println("# Found ", size(df_files_cleanest,1), " files considered clean 2.")

df_files_cleanest = df_files_use |>
    @filter( min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version ) |> 
    @filter( _.airmass <= args["max_airmass_clean"] ) |>
    @filter( abs(_.hour_angle) <= args["max_solar_hour_angle_clean"] ) |>
    @filter( _.expmeter_mean >= args["min_expmeter_clean"] ) |> 
    @filter( _.mean_pyroflux >= args["min_pyrhelio_clean"] ) |> 
    @filter( _.expmeter_rms <= args["max_expmeter_rms_frac_clean"]*_.expmeter_mean ) |> 
    @filter( _.rms_pyroflux <= args["max_pyrhelio_rms_frac_clean"]*_.mean_pyroflux  ) |> 
    @filter( _.expmeter_mean >= args["min_expmeter_to_pyrhelio_clean"]*_.mean_pyroflux ) |> 
    DataFrame
  println("# Found ", size(df_files_cleanest,1), " files considered clean 3.")

df_files_cleanest = df_files_use |>
    @filter( min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version ) |> 
    @filter( _.airmass <= args["max_airmass_clean"] ) |>
    @filter( abs(_.hour_angle) <= args["max_solar_hour_angle_clean"] ) |>
    @filter( _.expmeter_mean >= args["min_expmeter_clean"] ) |> 
    @filter( _.mean_pyroflux >= args["min_pyrhelio_clean"] ) |> 
    @filter( _.expmeter_rms <= args["max_expmeter_rms_frac_clean"]*_.expmeter_mean ) |> 
    @filter( _.rms_pyroflux <= args["max_pyrhelio_rms_frac_clean"]*_.mean_pyroflux  ) |> 
    @filter( _.expmeter_mean >= args["min_expmeter_to_pyrhelio_clean"]*_.mean_pyroflux ) |> 
    @filter(  extract_time_from_filename(_.Filename) >= start_time_clean ) |> 
    DataFrame
  println("# Found ", size(df_files_cleanest,1), " files considered clean 4.")

df_files_cleanest = df_files_use |>
    @filter( min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version ) |> 
    @filter( _.airmass <= args["max_airmass_clean"] ) |>
    @filter( abs(_.hour_angle) <= args["max_solar_hour_angle_clean"] ) |>
    @filter( _.expmeter_mean >= args["min_expmeter_clean"] ) |> 
    @filter( _.mean_pyroflux >= args["min_pyrhelio_clean"] ) |> 
    @filter( _.expmeter_rms <= args["max_expmeter_rms_frac_clean"]*_.expmeter_mean ) |> 
    @filter( _.rms_pyroflux <= args["max_pyrhelio_rms_frac_clean"]*_.mean_pyroflux  ) |> 
    @filter( _.expmeter_mean >= args["min_expmeter_to_pyrhelio_clean"]*_.mean_pyroflux ) |> 
    @filter(  extract_time_from_filename(_.Filename) <= stop_time_clean ) |> 
    DataFrame
  println("# Found ", size(df_files_cleanest,1), " files considered clean 5.")

df_files_cleanest = df_files_use |>
    @filter( min_drp_minor_version <= Base.thisminor(VersionNumber(_.drpversion)) <= max_drp_minor_version ) |> 
    @filter( _.airmass <= args["max_airmass_clean"] ) |>
    @filter( abs(_.hour_angle) <= args["max_solar_hour_angle_clean"] ) |>
    @filter( _.expmeter_mean >= args["min_expmeter_clean"] ) |> 
    @filter( _.mean_pyroflux >= args["min_pyrhelio_clean"] ) |> 
    @filter( _.expmeter_rms <= args["max_expmeter_rms_frac_clean"]*_.expmeter_mean ) |> 
    @filter( _.rms_pyroflux <= args["max_pyrhelio_rms_frac_clean"]*_.mean_pyroflux  ) |> 
    @filter( _.expmeter_mean >= args["min_expmeter_to_pyrhelio_clean"]*_.mean_pyroflux ) |> 
    @filter(  extract_time_from_filename(_.Filename) >= start_time_clean ) |> 
    @filter(  extract_time_from_filename(_.Filename) <= stop_time_clean ) |> 
    @orderby( abs(_.hour_angle) ) |> 
    @take(args["max_num_spectra_clean"] ) |>
    DataFrame
  println("# Found ", size(df_files_cleanest,1), " files considered clean.")


@assert size(df_files_cleanest,1) >= 1

  if hasproperty(df_files_cleanest,:dq1level)
    df_files_cleanest = df_files_cleanest |>   
       @filter( mod(_.dq1level,4) == 0  ) |>
       DataFrame
  end

if !(size(df_files_cleanest,1) >=1)
    @warn("No inputs passed all test for making clean spectra.")
end
@assert size(df_files_cleanest,1) >= 1

clean_obs_mask = map(fn->in(fn, df_files_cleanest.Filename),df_files_use.Filename) 

if verbose println(now()) end

if @isdefined(sed_filename)
   sed = Continuum.read_master_sed_neid(filename=sed_filename)
end
min_order_for_continuum = min_order(NEID2D()) # 12
max_order_for_continuum = max_order(NEID2D())-1  # Current DRP returns useless last order
orders_to_use_for_continuum = min_order_for_continuum:max_order_for_continuum

pipeline_plan = PipelinePlan()
 if verbose println("# Reading data files.")  end
 # For v1.0.0
 #all_spectra = Spectra2DBasic{Float64, Float32, Float32, Matrix{Float64}, Matrix{Float32}, Matrix{Float32}, NEID2D}[]
 # For v1.1.*
 all_spectra = Spectra2DBasic{Float64, Float64, Float32, Matrix{Float64}, Matrix{Float64}, Matrix{Float32}, NEID2D}[]
 spec = NEID.read_data(first(eachrow(df_files_use)).Filename, normalization=:blaze )
 mean_lambda = zeros(Float64,size(spec.flux));
 mean_clean_flux = zeros(Float64,size(spec.flux));
 mean_clean_var = zeros(Float64,size(spec.flux));
 mean_clean_flux_weight_sum = 0.0;
 mean_clean_flux_sed_normalized = zeros(Float64,size(spec.flux));
 mean_clean_var_sed_normalized = zeros(Float64,size(spec.flux));
 mean_clean_flux_sed_normalized_weight_sum = 0.0;
 mean_clean_flux_continuum_normalized = zeros(Float64,size(spec.flux));
 mean_clean_var_continuum_normalized = zeros(Float64,size(spec.flux));
 mean_clean_flux_continuum_normalized_weight_sum = 0.0;
 normalization_anchors_list = [];
 for (i,row) in enumerate(eachrow(df_files_use))
    if !(isfile(row.Filename)||islink(row.Filename))
        println("# Couldn't find input FITS file: ", row.Filename)
        continue
    end
    #= if !(isfile(row.continuum_filename)||islink(row.continuum_filename))
        println("# Couldn't find matching continuum file: ", row.continuum_filename)
        continue
    end =#
    if verbose
        println("# Reading ", row.Filename, "(",i,"/",size(df_files_use,1),")")
        flush(stdout)
    end
    local spec = NEID.read_data(row, normalization=:blaze )
    #file_hashes[row.Filename] = bytes2hex(sha256(row.Filename))
    file_hashes[row.Filename] = bytes2hex(open(md5,row.Filename))
    weight = 1
    if row.Filename ∈ df_files_cleanest.Filename
            mean_lambda .+= spec.λ
            mean_clean_flux .+= spec.flux # .*weight
            mean_clean_var .+= spec.var # .*weight
            global mean_clean_flux_weight_sum += weight
    end

    #=
    m = match(r"(neidL2_\d+[T_]\d+)\.fits$", row.Filename)
    continuum_filename = joinpath(neid_data_path,target_subdir,output_dir,"continuum", m[1] * "_continuum=new.jld2")
    if !(isfile(continuum_filename)||islink(continuum_filename))
        println("# Couldn't find matching continuum file", continuum_filename)
        continue
    end
    continuum = load(continuum_filename, "continuum")
    =#

    if @isdefined sed
      @warn("DRP v1.1 now provides blaze in L2 file.  This script has not been updated to use explicit an SED model.") 
      (f_norm, var_norm) = Continuum.normalize_by_sed(spec.λ,spec.flux,spec.var, sed; poly_order = args["order_poly_continuum"], half_width = args["continuum_poly_half_width"], quantile = args["quantile_fit_continuum"], orders_to_use=orders_to_use_for_continuum)
      sed_norm_failed = any([all(isnan.(view(f_norm,:,ord))) for ord in 1:(size(f_norm,2)-1)]) #check if any orders besides the last order are all nans
      if row.Filename ∈ df_files_cleanest.Filename && !sed_norm_failed
            mean_clean_flux_sed_normalized .+= f_norm # .*weight
            mean_clean_var_sed_normalized .+= var_norm # .*weight
            global mean_clean_flux_sed_normalized_weight_sum +=  weight
      end
    else
        f_norm = spec.flux
        var_norm = spec.var
    end

 
   if args["apply_continuum_normalization"] && args["continuum_normalization_individually"]
        println("Applying continuum normalization to each spectrum individually.")
        local anchors, continuum, f_filtered
        if args["anchors_filename"] != nothing
            @assert isfile(args["anchors_filename"]) && filesize(args["anchors_filename"])>0
            anchors = load(args["anchors_filename"],"anchors")
            (anchors, continuum, f_filtered) = Continuum.calc_continuum(spec.λ, f_norm, var_norm, anchors;
                orders_to_use = orders_to_use_for_continuum, verbose = true )
            push!(normalization_anchors_list, anchors)
        else
            (anchors, continuum, f_filtered) = Continuum.calc_continuum(spec.λ, f_norm, var_norm; fwhm = args["fwhm_continuum"]*1000, ν = args["nu_continuum"],
                stretch_factor = args["stretch_factor"], merging_threshold = args["merging_threshold"], smoothing_half_width = args["smoothing_half_width"], min_R_factor = args["min_rollingpin_r"],
                orders_to_use = orders_to_use_for_continuum, verbose = true )
            push!(normalization_anchors_list, anchors)
        end
            #=
            continuum = load(row.continuum_filename, "continuum")
            #file_hashes[row.continuum_filename] = bytes2hex(sha256(row.continuum_filename))
            file_hashes[row.continuum_filename] = bytes2hex(open(md5,row.continuum_filename))
            =#
        #f_norm ./= continuum
        #var_norm ./= continuum.^2
        for ord in orders_to_use_for_continuum
	    f_norm[:,ord] ./= view(continuum,:,ord)
            var_norm[:,ord] ./= view(continuum,:,ord).^2
	end
        #continuum_norm_failed = any([all(isnan.(all_spectra[i].flux[:,i])) for i in 1:size(all_spectra[i].flux,2)][1:end-1]) #check if any orders besides the last order are all nans
        continuum_norm_failed = any([all(isnan.(view(f_norm,:,ord))) for ord in orders_to_use_for_continuum]) #check if any orders used for continuum are all nans
        if row.Filename ∈ df_files_cleanest.Filename && !continuum_norm_failed
                mean_clean_flux_continuum_normalized .+= f_norm # .*weight
                mean_clean_var_continuum_normalized .+= var_norm # .*weight
                global mean_clean_flux_continuum_normalized_weight_sum += weight
        end
        spec.flux .= f_norm
        spec.var .= var_norm

    end
    push!(all_spectra,spec)
 end
 GC.gc()
 mean_lambda ./= mean_clean_flux_weight_sum
 mean_clean_flux ./= mean_clean_flux_weight_sum
 mean_clean_var ./= mean_clean_flux_weight_sum
 mean_clean_flux_sed_normalized ./= mean_clean_flux_sed_normalized_weight_sum
 mean_clean_var_sed_normalized ./= mean_clean_flux_sed_normalized_weight_sum
 dont_need_to!(pipeline_plan,:read_spectra);


 if args["apply_continuum_normalization"] && !args["continuum_normalization_individually"]
     println("Applying continuum normalization based on mean of clean spectra.")
     local anchors, continuum, f_filtered
     if args["anchors_filename"] !=nothing
         @assert isfile(args["anchors_filename"]) && filesize(args["anchors_filename"])>0
         anchors = load(args["anchors_filename"],"anchors")
     else
         println("# Computing Continuum model.")
         if @isdefined sed
             @warn("DRP v1.1 now provides blaze in L2 file.  This script has not been updated to use explicit an SED model.") 
             (anchors, continuum, f_filtered) = Continuum.calc_continuum(spec.λ, mean_clean_flux_sed_normalized, mean_clean_var_sed_normalized; fwhm = args["fwhm_continuum"]*1000, ν = args["nu_continuum"],
                stretch_factor = args["stretch_factor"], merging_threshold = args["merging_threshold"], smoothing_half_width = args["smoothing_half_width"], min_R_factor = args["min_rollingpin_r"],
                orders_to_use = orders_to_use_for_continuum, verbose = false )
            if !isnothing(args["anchors_filename_output"])
                println("# Storing anchors used for continuum model in ",args["anchors_filename_output"],  ".")
                save(args["anchors_filename_output"], Dict("anchors" => anchors) )
            end
        else
            (anchors, continuum, f_filtered) = Continuum.calc_continuum(spec.λ, mean_clean_flux, mean_clean_var; fwhm = args["fwhm_continuum"]*1000, ν = args["nu_continuum"],
                stretch_factor = args["stretch_factor"], merging_threshold = args["merging_threshold"], smoothing_half_width = args["smoothing_half_width"], min_R_factor = args["min_rollingpin_r"],
                orders_to_use = orders_to_use_for_continuum, verbose = true )
            if !isnothing(args["anchors_filename_output"])
                println("# Storing anchors used for continuum model in ",args["anchors_filename_output"],  ".")
                save(args["anchors_filename_output"], Dict("anchors" => anchors) )
            end
        end # @isdefined sed
    end # args["anchors_filename"] 
    normalization_anchors_list = anchors

    weight = 1
    for (i,row) in enumerate(eachrow(df_files_use))
        (anchors, continuum, f_filtered) = Continuum.calc_continuum(all_spectra[i].λ, all_spectra[i].flux, all_spectra[i].var,
                    anchors, smoothing_half_width = args["smoothing_half_width"], orders_to_use=orders_to_use_for_continuum)
        for ord in orders_to_use_for_continuum
	    all_spectra[i].flux[:,ord] ./= view(continuum,:,ord)
            all_spectra[i].var[:,ord] ./= view(continuum,:,ord).^2
	end
        continuum_norm_failed = any([all(isnan.(view(all_spectra[i].flux,:,ord))) for ord in orders_to_use_for_continuum]) #check if any orders used for continuum are all nans
        if row.Filename ∈ df_files_cleanest.Filename && !continuum_norm_failed
            mean_clean_flux_continuum_normalized .+= all_spectra[i].flux # .*weight
            mean_clean_var_continuum_normalized .+= all_spectra[i].var # .*weight
            global mean_clean_flux_continuum_normalized_weight_sum += weight
        end
    end
 end
 mean_clean_flux_continuum_normalized ./= mean_clean_flux_continuum_normalized_weight_sum
 mean_clean_var_continuum_normalized ./= mean_clean_flux_continuum_normalized_weight_sum

 order_list_timeseries = extract_orders(all_spectra, pipeline_plan,  orders_to_use=orders_to_use, remove_bad_chunks=false, recalc=true )


line_width = line_width_50_default
 max_mask_scale_factor = max(mask_scale_factor, max(lsf_width,
                             line_width/sqrt(8*log(2)))/default_ccf_mask_v_width(NEID2D())) # 4.0
 max_bc = RvSpectMLBase.max_bc
 if args["target"] == Sun || args["target"] == "Solar"
    max_bc = RvSpectMLBase.max_bc_solar
 end
 #max_orders = min_order(NEID2D()):max_order(NEID2D())
 #good_orders = orders_to_use_default(NEID2D())
 #orders_to_use = max_orders
 if isfile(line_list_filename) && (filesize(line_list_filename)>0) && !args["recompute_line_weights"]
   println("# Reading ", line_list_filename)
   line_list_espresso = CSV.read(line_list_filename, DataFrame)
   @assert all(map(k->k ∈ names(line_list_espresso), ["lambda","weight","order"]))
   dont_need_to!(pipeline_plan,:clean_line_list_tellurics)
 else
    #orders_to_use = good_orders
    #order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=orders_to_use, recalc=true )
    touch(line_list_filename)
    if (isfile(line_list_filename) && (filesize(line_list_filename)>0))
        line_list_input_filename = line_list_filename
    else
        line_list_input_filename = joinpath(pkgdir(EchelleCCFs),"data","masks","espresso+neid_mask_97_to_108.mas")
    end
    println("# Recreating line list weights from ", line_list_input_filename)
    line_list_espresso = prepare_line_list(line_list_input_filename, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity,
       Δv_to_avoid_tellurics = 2*max_bc+range_no_mask_change*line_width_50_default+max_mask_scale_factor*default_ccf_mask_v_width(NEID2D()), orders_to_use=orders_to_use, recalc=true, verbose=true)
     if args["recompute_line_weights"] && !isnothing(args["line_list_output_filename"])
        CSV.write(args["line_list_output_filename"], line_list_espresso)
     end
 end
 #file_hashes[line_list_filename] = bytes2hex(sha256(line_list_filename))
 file_hashes[line_list_filename] = bytes2hex(open(md5,line_list_filename))
 #outputs["line_list_espresso"] = line_list_espresso

if args["variable_mask_scale"]
     #maxΔfwhm² = -0.569375
     maxΔfwhm² = -0.56930
     @assert maximum(df_files_use.Δfwhm²) <= maxΔfwhm²
     if maximum(df_files_use.Δfwhm²) > maxΔfwhm²
        println("# Warning: dangerously large maximum(df_files_use.Δfwhm²) = ", maximum(df_files_use.Δfwhm²a) )
     end
     Δfwhm = 1000.0 .* sqrt.(clamp.(maxΔfwhm².-df_files_use.Δfwhm²[1:length(all_spectra)], 0., Inf))  # How much to increase fwhm by to acheive uniform fwhm
else
    Δfwhm = zeros(0)
end

if verbose println(now()) end
line_width_50 = line_width_50_default
  #order_list_timeseries = extract_orders(all_spectra,pipeline_plan,  orders_to_use=orders_to_use,  remove_bad_chunks=false, recalc=true )
  @time (order_ccfs, order_ccf_vars, v_grid_order_ccfs) = ccf_orders(order_list_timeseries, line_list_espresso, pipeline_plan,
    mask_type=:gaussian, Δfwhm=Δfwhm,
    mask_scale_factor=mask_scale_factor, range_no_mask_change=5*line_width_50, ccf_mid_velocity=ccf_mid_velocity, v_step=v_step,
    v_max=max(range_no_mask_change*line_width_50,2*max_bc), orders_to_use=orders_to_use, allow_nans=true, calc_ccf_var=true,
    recalc=true)

#orders_to_use2 = orders_to_use[map(i->!iszero(order_ccfs[:,i,:]),1:length(orders_to_use))]
#order_list_timeseries2 = extract_orders(all_spectra,pipeline_plan,  orders_to_use=orders_to_use2,  remove_bad_chunks=false, recalc=true )

#=
ccf_dir = joinpath(neid_data_path,target_subdir,"ccfs")
if !isdir(ccf_dir)
   mkdir(ccf_dir)
end
=#
if verbose println(now()) end
#daily_ccf_filename = joinpath(neid_data_path,target_subdir,"output","ccfs", date_str * "_ccfs=default.jld2")
#daily_ccf_filename = joinpath(neid_data_path,target_subdir,"ccfs", date_str * "_ccfs=default.jld2")
println("# Saving results to ", daily_ccf_filename, ".")
  stop_processing_time = now()
  jldopen(daily_ccf_filename, "w") do f
    f["v_grid"] = collect(v_grid_order_ccfs)
    f["order_ccfs"] = order_ccfs
    f["order_ccf_vars"] = order_ccf_vars
    f["Δfwhm"] = Δfwhm 
    f["orders_to_use"] = orders_to_use
    f["manifest"] = df_files_use
    f["clean_obs_mask"] = clean_obs_mask
    f["calc_order_ccf_args"] = args
    f["ccf_line_list"] = line_list_espresso
    if (size(df_files_cleanest,1) >= 1) && any(.!iszero.(mean_clean_flux))
      f["mean_lambda"] = mean_lambda
      f["mean_clean_flux"] = mean_clean_flux
      f["mean_clean_var"] = mean_clean_var
      if (size(mean_clean_flux_sed_normalized,1) > 0) && any(.!iszero.(mean_clean_flux_sed_normalized)) && any(.!isnan.(mean_clean_flux_sed_normalized))
          f["mean_clean_flux_sed_normalized"] = mean_clean_flux_sed_normalized
          f["mean_clean_var_sed_normalized"] = mean_clean_var_sed_normalized
      end
      if (size(mean_clean_flux_continuum_normalized,1) > 0) && any(.!iszero(mean_clean_flux_continuum_normalized))  && any(.!isnan.(mean_clean_flux_continuum_normalized))
          f["mean_clean_flux_continuum_normalized"] = mean_clean_flux_continuum_normalized
          f["mean_clean_var_continuum_normalized"] = mean_clean_var_continuum_normalized
      end
    end
    if args["apply_continuum_normalization"]
        f["normalization_anchors"] = normalization_anchors_list
    end
    f["start_processing_time"] = start_processing_time
    f["stop_processing_time"] = stop_processing_time
    f["daily_ccf_filename"] = daily_ccf_filename
    f["drpversion"] = max_drp_minor_version 
    f["file_hashes"] = file_hashes
  end
if verbose println(now()) end

#=
for (i,row) in enumerate(eachrow(df_files_use))
    m = match(r"(neidL2_\d+[T_]\d+)\.fits$", row.Filename)
    #ccf_filename = joinpath(neid_data_path,target_subdir,"output","ccfs", m[1] * "_ccfs=default.jld2")
    ccf_filename = joinpath(neid_data_path,target_subdir,"ccfs", m[1] * "_ccfs=default.jld2")
    jldopen(ccf_filename, "w") do f
        f["v_grid"] = collect(v_grid_order_ccfs)
        f["order_ccfs"] = order_ccfs[:,:,i]
        f["order_ccf_vars"] = order_ccf_vars[:,:,i]
        f["orders_to_use"] = orders_to_use
        for k in keys(row)
            f[string(k)] = row[k]
        end
        f["ccf_filename"] = ccf_filename
    end
    if i==size(order_ccfs,3) break end
end
=#


#=

#good_orders = orders_to_use_default(NEID2D())
 #good_orders = 56:111 # orders_to_use_default(NEID2D())
 #good_orders = 57:96 # orders_to_use_default(NEID2D())
 #good_orders = 55 .+ vcat(1:14, 17:24, 27:29, 37:38,40:41)
 #good_orders = 55 .+ collect(1:14)
 #good_orders = 55 .+ collect(17:24)
 #good_orders = 55 .+ vcat(collect(17:24), collect(27:29), collect(37:38),collect(40:41))
 #good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41))
 #good_orders = [55, 56, 57]
# order_list_timeseries = extract_orders(all_spectra,pipeline_plan,  orders_to_use=good_orders, remove_bad_chunks=false, recalc=true )

outputs["times"] = order_list_timeseries.times
 alt_sun = map(i->calc_solar_alt(all_spectra[i].metadata[:bjd]),1:length(all_spectra))
 obs_idx_max_alt_sun = argmax(alt_sun)
 time_max_alt_sun = order_list_timeseries.times[obs_idx_max_alt_sun]
 outputs["time_max_alt_sun"] = time_max_alt_sun

 order_snrs = [RvSpectMLBase.calc_snr(order_list_timeseries[obsid].data[c]) for c in 1:length(order_list_timeseries[1].data),
                    obsid in 1:length(order_list_timeseries) ]
 order_sum_flux = [NaNMath.sum(order_list_timeseries[obsid].data[c].flux) for c in 1:length(order_list_timeseries[1].data),
                                        obsid in 1:length(order_list_timeseries) ]
 order_mean_flux = [NaNMath.mean(order_list_timeseries[obsid].data[c].flux) for c in 1:length(order_list_timeseries[1].data),
                                                                        obsid in 1:length(order_list_timeseries) ]

 mean_order_mean_flux = vec(mean(order_mean_flux,dims=2))


 order_weights =  (order_snrs[:,obs_idx_max_alt_sun]./order_snrs).^2
 #order_weights = (order_sum_flux[:,obs_idx_max_alt_sun]./order_sum_flux)
 outputs["order_snrs"] = order_snrs
 outputs["order_sum_flux"] = order_sum_flux
 outputs["order_mean_flux"] = order_mean_flux
 outputs["order_weights"] = order_weights

#times = order_list_timeseries.times.-time_max_alt_sun
#plot(times,order_weights', label=:none)
=#

#=
for obsid in 1:length(order_list_timeseries)
    for ch in 1:length(order_list_timeseries[1].data)
        order_list_timeseries[obsid].data[ch].flux .*= order_weights[ch,obsid] #./mean_order_mean_flux[ch]
        order_list_timeseries[obsid].data[ch].var .*= (order_weights[ch,obsid]).^2 # ./mean_order_mean_flux[ch]).^2
    end
end
=#

#=for obsid in 1:length(order_list_timeseries)
    for ch in 1:length(order_list_timeseries[1].data)
        order_list_timeseries[obsid].data[ch].flux ./= order_weights[ch,obsid] #./mean_order_mean_flux[ch]
        order_list_timeseries[obsid].data[ch].var ./= (order_weights[ch,obsid]).^2 # ./mean_order_mean_flux[ch]).^2
    end
end
=#

#=
using Plots
 plot(map(obsid->NaNMath.sum(order_list_timeseries[obsid].data[1].flux)/(NaNMath.sum(order_list_timeseries[obsid].data[end].flux)),
                    1:length(all_spectra)),label="Raw")
 #plot!(map(obsid->NaNMath.sum(order_list_timeseries[obsid].data[10].flux)*order_weights[10,obsid]/(NaNMath.sum(order_list_timeseries[obsid].data[30].flux*order_weights[30,obsid])), 1:length(all_spectra)),label="weighted")
 #plot!(map(obsid->NaNMath.sum(order_list_timeseries[obsid].data[10].flux)*order_weights[10,obsid]^2/(NaNMath.sum(order_list_timeseries[obsid].data[30].flux*order_weights[30,obsid]^2)), 1:length(all_spectra)),label="weighted^2")

# plot(order_list_timeseries[obsid].data[ch].λ,order_list_timeseries[obsid].data[ch].flux)
=#

#=
times = order_list_timeseries.times.-time_max_alt_sun
for obsid in 1:length(order_list_timeseries)
    order_snrs[:,obsid] ./= order_snrs[15,obsid]
    order_snrs_new[:,obsid] ./= order_snrs_new[15,obsid]
end

plot(times,order_snrs', label=:none)
plot!(times,order_snrs_new'./order_snrs_new[15,:]', label=:none)

plt = plot()
for obsid in 1:length(order_list_timeseries)
    for ch in 1:10:length(order_list_timeseries[1].data)
        plot(order_list_timeseries[obsid].data[ch].λ,order_list_timeseries[obsid].data[ch].flux)

    end
end
=#

#=
maxΔfwhm² = -0.569375
@assert maximum(df_files_use.Δfwhm²) < maxΔfwhm²
Δfwhm = 1000.0 .* sqrt.(maxΔfwhm².-df_files_use.Δfwhm²)  # How much to increase fwhm by to acheive uniform fwhm

#msf = 1.4; fwtf = 1.5
#msf = 6.0; fwtf = 1.5
#msf = 0.62; fwtf = 1.5  # one pixel because tophat not using NEID's default v width
#msf = 1.0; fwtf = 1.5  # using NEID's default v width
msf = lsf_width/default_ccf_mask_v_width(NEID2D()); fwtf = 0.5  # using LSF width
#msf = line_width_50_default/default_ccf_mask_v_width(NEID2D()); fwtf = 1.0  # using line width
#msf = 1.24; fwtf = 1.5   # two pixels
 println("# mask_scale_factor = ", msf, ", frac_width_to_fit =  ", fwtf, ". ")
 ((ccfs_espresso, ccf_vars_espresso), v_grid) = ccf_total(order_list_timeseries, line_list_espresso, pipeline_plan,
   mask_type=:gaussian, Δfwhm=Δfwhm,
   mask_scale_factor=msf, range_no_mask_change=5*line_width_50, ccf_mid_velocity=ccf_mid_velocity, v_step=100, #155,
   v_max=max(5*line_width_50,2*max_bc), allow_nans=true, calc_ccf_var=true, recalc=true)
  outputs["v_grid"] = collect(v_grid)
  outputs["ccfs_espresso"] = ccfs_espresso
  outputs["ccf_vars_espresso"] = ccf_vars_espresso

ccf_median = vec(median(ccfs_espresso,dims=2))
 ccf_dist_from_median = map(i->RvSpectML.earth_mover_distance(v_grid,ccf_median,ccfs_espresso[:,i]), 1:size(ccfs_espresso,2) )

fwtf = 0.5; println("# mask_scale_factor = ", msf, ", frac_width_to_fit =  ", fwtf, ". ")
 alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
 rvs_ccf_espresso = RvSpectML.calc_rvs_from_ccf_total(ccfs_espresso, ccf_vars_espresso, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
 σ_rvs_espresso = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))
 #scatter((order_list_timeseries.times.-first(order_list_timeseries.times)).*24, rvs_ccf_espresso, label=:none)
 outputs["rvs_ccf_espresso"] = rvs_ccf_espresso
 outputs["σ_rvs_espresso"] = σ_rvs_espresso

 delay_until_first_obs_used_for_fit = 1/24
 times = copy(order_list_timeseries.times)
 #Δt_fit = min(time_max_alt_sun-(first(times)+delay_until_first_obs_used_for_fit),last(times)-time_max_alt_sun)
 Δt_fit = 2.0/24
 idx_first_obs_to_fit = findfirst(t->time_max_alt_sun-t<=Δt_fit,times)
 @assert 1 <= idx_first_obs_to_fit length(times)
 idx_last_obs_to_fit = findlast(t->t-time_max_alt_sun<=Δt_fit,times)
 idx_to_fit = idx_first_obs_to_fit:idx_last_obs_to_fit
 idx_to_fit = idx_to_fit[ccf_dist_from_median[idx_to_fit] .< 10]
 (mean_rv_to_fit, std_rv_to_fit) = mean_and_std(rvs_ccf_espresso[idx_to_fit])
 idx_for_robust_fit = idx_to_fit[mean_rv_to_fit-3*std_rv_to_fit .<= rvs_ccf_espresso[idx_to_fit] .<= mean_rv_to_fit+3*std_rv_to_fit]

 times .-= time_max_alt_sun
 #times = order_list_timeseries.times.-first(order_list_timeseries.times)
 x = [times[idx_for_robust_fit]  ones(length(times[idx_for_robust_fit])) ]
 linfit = (x'*x)\x'*(rvs_ccf_espresso[idx_for_robust_fit])
 x = [times  ones(length(times)) ]
 v_fit = (x*linfit)
 v_t0 = [[0] [1]]*linfit
 idx_not_bad_distance = ccf_dist_from_median .< 10
 (times_binned, rvs_binned) = bin_times_and_rvs_max_Δt(times=times[idx_not_bad_distance], rvs=rvs_ccf_espresso[idx_not_bad_distance], Δt_threshold=5/(60*24))
 #idx_to_fit_binned = findfirst(t->t>=times[first(idx_to_fit)], times_binned):findlast(t->t<=times[last(idx_to_fit)], times_binned)
 #idx_to_fit_binned = findfirst(t->time_max_alt_sun-t.<=Δt_fit,times_binned):findlast(t->t-time_max_alt_sun.<=Δt_fit,times_binned)
 idx_to_fit_binned = findfirst(t->t>=-Δt_fit,times_binned):findlast(t->t<=Δt_fit,times_binned)
 x = [times_binned  ones(length(times_binned)) ]
 v_fit_binned = x*linfit
 println("# RMS to linear fit: ", std(rvs_ccf_espresso[idx_to_fit].-v_fit[idx_to_fit]), " m/s.  Binned: ", std(rvs_binned[idx_to_fit_binned].-v_fit_binned[idx_to_fit_binned]))
 println("# Mean RV = ", mean(rvs_binned[idx_to_fit_binned]), " m/s  Slope = ", linfit[1]/24, " m/s/hr")
 outputs["times_binned"] = times_binned
 outputs["rvs_binned"] = rvs_binned
 outputs["rv_rms"] = std(rvs_ccf_espresso[idx_to_fit].-v_fit[idx_to_fit])
 outputs["rv_0"] = v_t0
 outputs["rv_mean"] = mean(rvs_binned[idx_to_fit_binned])
 outputs["rv_slope"] = linfit[1]/24
 outputs["rv_rms_binned"] = std(rvs_binned[idx_to_fit_binned].-v_fit_binned[idx_to_fit_binned])

using Plots
 plt = scatter( (times[idx_not_bad_distance]).*24, (rvs_ccf_espresso.-v_fit)[idx_not_bad_distance],  ms=1.5, label=:none)
 scatter!(plt,(times_binned).*24, (rvs_binned.-v_fit_binned), ms=4.0, label=:none)
 vline!(plt,[-Δt_fit,Δt_fit].* 24, label=:none, color=:grey)
 vline!(plt,[0].* 24, label=:none, color=:orange)
 xlims!(plt,-3.5,3.5)
 #ylims!(plt,minimum(rvs_ccf_espresso.-v_fit)-0.5*std_rv_to_fit, maximum(rvs_ccf_espresso.-v_fit)+0.5*std_rv_to_fit)
 ylims!(plt, extrema((rvs_ccf_espresso.-v_fit)[idx_not_bad_distance]))
 xlabel!(plt,"Time (hr)")
 ylabel!(plt,"RV - linear fit (m/s)")
 title!(plt,"Sun " * date_str * " CCF ESPRESSO\n" *
        "rms_5m = " * string(round(std(rvs_binned[idx_to_fit_binned].-v_fit_binned[idx_to_fit_binned]),digits=3)) * "m/s" *
        " mean = " * string(round(mean(rvs_binned[idx_to_fit_binned]),digits=2)) *
        " Slope = " * string(round(linfit[1]/24,digits=2))  * "m/s/hr")
 display(plt)
 savefig(joinpath(output_dir,"rvs_" * date_str * ".png"))


(ccfs_norm, ccf_vars_norm) = EchelleCCFs.calc_normalized_ccfs( v_grid, ccfs_espresso, ccf_vars_espresso, v_mid=ccf_mid_velocity, dv_min=2.5*line_width_50,dv_max= 2*max_bc )
 ccf_template = EchelleCCFs.calc_ccf_template(ccfs_norm[:,idx_to_fit], ccf_vars_norm[:,idx_to_fit], assume_normalized=true)
 ccf_vars_template = EchelleCCFs.calc_ccf_template(ccf_vars_norm[:,idx_to_fit], ccf_vars_norm[:,idx_to_fit], assume_normalized=true)
 #plot(v_grid,ccfs_espresso, labels=:none)
 #plot(v_grid,ccfs_norm , labels=:none)
 #plt = plot(v_grid,(ccfs_norm.-ccf_template)[:,range(1,size(ccfs_norm,2),step=2)], labels=:none)
 #xlims!(plt,ccf_mid_velocity-2*line_width_50,ccf_mid_velocity+2*line_width_50)
 #heatmap(v_grid,1:size(ccfs_norm,2),(ccfs_norm.-ccf_template)',seriescolor=:balance)
 cols_to_fit = findlast(x->x<ccf_mid_velocity-3*line_width_50,v_grid):findfirst(x->x>ccf_mid_velocity+3*line_width_50,v_grid)
 maxabsz = maximum(abs.(extrema(view(ccfs_norm,cols_to_fit,idx_to_fit)'.-view(ccf_template,cols_to_fit,:)')))
 plt = heatmap(view(v_grid,cols_to_fit),(times[1:size(ccfs_norm,2)]).*24,view(ccfs_norm,cols_to_fit,:)'.-view(ccf_template,cols_to_fit,:)',clims=(-maxabsz,maxabsz),seriescolor=:balance)
 xlabel!(plt,"v (m/s)")
 ylabel!(plt,"Observation Time (hr)")
 title!(plt,"Sun " * date_str * ":\nCCF - <CCF>, ESPRESSO mask + SNR weights")
 #xlims!(ccf_mid_velocity-RvSpectMLBase.max_bc,ccf_mid_velocity+RvSpectMLBase.max_bc)
 xlims!(plt,ccf_mid_velocity-2*line_width_50,ccf_mid_velocity+2*line_width_50)
 hline!(plt,[-Δt_fit,Δt_fit].*24, label=:none, color=:grey)
 hline!(plt,[0.], label=:none, color=:orange)
 ylims!(plt,-3.5,3.0)
 #vline!([-2.1e4], color=4, label=:none)
 display(plt)
 #savefig(joinpath(output_dir,"ccf_heatmaps_sun_" * date_str * "_espresso_msf=" * string(round(msf,digits=1)) * ".png"))
 savefig(joinpath(output_dir,"ccf_heatmaps_sun_" * date_str * "_espresso_msf=" * string(round(msf,digits=1)) * "_new.png"))

 outputs["mean_ccf"] = ccf_template

gp_post = RvSpectML.TemporalGPInterpolation.construct_gp_posterior(v_grid,ccf_template,sigmasq_obs=ccf_vars_template,use_logx=false,smooth_factor=2)
 ccf_template_smooth = RvSpectML.TemporalGPInterpolation.predict_mean(gp_post(v_grid))
 ccf_template_deriv = RvSpectML.numerical_deriv(v_grid,RvSpectML.TemporalGPInterpolation.predict_mean(gp_post(v_grid)))
 norm = sum(abs2.(ccf_template_deriv[cols_to_fit]))
 rvs_proj = vec(sum((ccf_template_smooth[cols_to_fit].-ccfs_norm[cols_to_fit,:]).*ccf_template_deriv[cols_to_fit],dims=1)./norm)
 σ_rvs_proj = vec(sum((ccf_vars_norm[cols_to_fit]).*abs2.(ccf_template_deriv[cols_to_fit]),dims=1))./ norm
 #ccf_resid_minus_rv_proj = (ccfs_norm.-ccf_template).+(rvs_proj.*ccf_template_deriv')'
 ccf_resid_minus_rv_proj = (ccfs_norm.-ccf_template).+((v_fit.-mean(v_fit)).*ccf_template_deriv')'

 outputs["ccf_template_smooth"] = ccf_template_smooth
 outputs["rvs_proj"] = rvs_proj
 outputs["σ_rvs_proj"] = σ_rvs_proj
 outputs["ccf_resid_minus_rv_proj"] = ccf_resid_minus_rv_proj

#plt = heatmap(view(v_grid,cols_to_fit),(times[1:size(ccfs_norm,2)]).*24,ccf_resid_minus_rv_proj[cols_to_fit,:]')
maxabsz = maximum(abs.(extrema(ccf_resid_minus_rv_proj)))
climmax = 0.00030
 plt = heatmap(v_grid,(times[1:size(ccfs_norm,2)]).*24,ccf_resid_minus_rv_proj',clims=(-climmax,climmax),seriescolor=:balance)
 xlabel!(plt,"v (m/s)")
 ylabel!(plt,"Observation Time (hr)")
 title!(plt,"Sun " * date_str * ":\nCCF-<CCF>+v⋅d<CCF>/dv")
 #xlims!(ccf_mid_velocity-RvSpectMLBase.max_bc,ccf_mid_velocity+RvSpectMLBase.max_bc)
 xlims!(plt,ccf_mid_velocity-2*line_width_50,ccf_mid_velocity+2*line_width_50)
 hline!(plt,[-Δt_fit,Δt_fit].*24, label=:none, color=:grey)
 hline!(plt,[0.], label=:none, color=:orange)
 ylims!(plt,-3.5,3.0)
 display(plt)
 #savefig(joinpath(output_dir,"ccf_resid_heatmaps_sun_" * date_str * "_espresso_msf=" * string(round(msf,digits=1)) * ".png"))
 savefig(joinpath(output_dir,"ccf_resid_heatmaps_sun_" * date_str * "_espresso_msf=" * string(round(msf,digits=1)) * "_new.png"))
=#
@assert true

#=
 all_spectra = nothing
 order_list_timeseries = nothing
 GC.gc()
end
=#

#EchelleCCFs.write_each_ccf_fits(map(s->s.metadata,all_spectra), v_grid, ccfs_espresso, ccf_vars_espresso,   fits_hdr=EchelleCCFs.make_ccf_fits_header(all_spectra[1].metadata),
#            line_list_filename=line_list_filename,total_vel=rvs_ccf_espresso, sigma_total_vel=σ_rvs_espresso )

#=
##
max_orders = 56:111
  order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=max_orders, recalc=true )

line_list_filename = "G2.espresso.mas"
 linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",line_list_filename)
 line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity,
      Δv_to_avoid_tellurics = 1*RvSpectMLBase.max_bc+5*line_width_50+2*default_ccf_mask_v_width(EXPRES1D()),
      orders_to_use=max_orders, recalc=true, verbose=true)

##

good_orders = orders_to_use_default(NEID2D())

(order_ccfs, order_ccf_vars, v_grid_order_ccfs) = ccf_orders(order_list_timeseries, line_list_espresso, pipeline_plan, mask_scale_factor=msf, ccf_mid_velocity=ccf_mid_velocity,
 v_step=155, v_max=max(5*line_width_50,2*max_bc), orders_to_use=good_orders, allow_nans=true, calc_ccf_var=true,recalc=true)
 outputs["order_ccfs"] = order_ccfs
 outputs["order_ccf_vars"] = order_ccf_vars
 outputs["v_grid_order_ccfs"] =v_grid_order_ccfs

obs_id = idx_to_fit[1]
 plt = heatmap(v_grid_order_ccfs, 1:size(order_ccfs,2), (order_ccfs[:,:,obs_id]./mean(order_ccfs[:,:,obs_id],dims=1))')
 xlabel!(plt,"v (m/s)")
 ylabel!(plt,"Order (ignore #)")
 title!(plt,"Sun " * date_str  * " Obs# " * string(obs_id))
 #xlims!(ccf_mid_velocity-RvSpectMLBase.max_bc,ccf_mid_velocity+RvSpectMLBase.max_bc)
 xlims!(plt,ccf_mid_velocity-2.0*line_width_50,ccf_mid_velocity+2.0*line_width_50)
 display(plt)
 savefig(joinpath(output_dir,"order_ccf_heatmaps_sun_" * date_str * "_espresso_msf=" * string(round(msf,digits=1)) * ".png"))


ccf_order_template = zeros(size(order_ccfs,1),size(order_ccfs,2))
 ccf_order_template_vars = zeros(size(order_ccfs,1),size(order_ccfs,2))
 ccf_order_resid = zeros(size(order_ccfs,1),size(order_ccfs,2),size(order_ccfs,3))
 rvs_proj_order = zeros(length(times),size(order_ccfs,2))
 σ_rvs_proj_order = zeros(length(times),size(order_ccfs,2))
 ccf_order_resid_rv_from_order = zeros(size(order_ccfs,1),size(order_ccfs,2),size(order_ccfs,3))

 for ord in 1:size(order_ccfs,2)
     if !(sum(order_ccfs[:,ord,:]) > 0)  continue end
     (ccfs_order_norm, ccf_order_vars_norm) = EchelleCCFs.calc_normalized_ccfs( v_grid, order_ccfs[:,ord,:], order_ccf_vars[:,ord,:], v_mid=ccf_mid_velocity, dv_min=2.5*line_width_50,dv_max= 2*max_bc )
     ccf_order_template[:,ord] .= EchelleCCFs.calc_ccf_template(ccfs_order_norm[:,idx_to_fit], ccf_order_vars_norm[:,idx_to_fit], assume_normalized=true)
     ccf_order_template_vars[:,ord] .= EchelleCCFs.calc_ccf_template(ccf_order_vars_norm[:,idx_to_fit], ccf_order_vars_norm[:,idx_to_fit], assume_normalized=true)
     println("# ord = ", ord, " mean = ", NaNMath.mean(ccf_order_template[:,ord]), " var = ", NaNMath.mean(ccf_order_template_vars[:,ord]))
     gp_post = RvSpectML.TemporalGPInterpolation.construct_gp_posterior(v_grid,ccf_order_template[:,ord],sigmasq_obs=ccf_order_template_vars[:,ord],use_logx=false,smooth_factor=2)
     ccf_template_smooth = RvSpectML.TemporalGPInterpolation.predict_mean(gp_post(v_grid))
     ccf_template_deriv = RvSpectML.numerical_deriv(v_grid,RvSpectML.TemporalGPInterpolation.predict_mean(gp_post(v_grid)))
     ccf_order_resid[:,ord,:] .= (ccfs_order_norm.-ccf_order_template[:,ord]).+(rvs_proj.*ccf_template_deriv')'
     norm = sum(abs2.(ccf_template_deriv[cols_to_fit]))
     rvs_proj_order[:,ord] .= vec(sum((ccf_template_smooth[cols_to_fit].-ccfs_order_norm[cols_to_fit,:]).*ccf_template_deriv[cols_to_fit],dims=1))./norm
     σ_rvs_proj_order[:,ord] .= vec(sum((ccf_order_vars_norm[cols_to_fit,:]).*abs2.(ccf_template_deriv[cols_to_fit]),dims=1))./ norm
     ccf_order_resid_rv_from_order[:,ord,:] .= (ccfs_order_norm.-ccf_template).+(rvs_proj_order[:,ord].*ccf_template_deriv')'
 end
 outputs["ccf_order_template"] = ccf_order_template
 outputs["ccf_order_template_vars"] = ccf_order_template_vars
 outputs["rvs_proj_order"] = rvs_proj_order
 outputs["σ_rvs_proj_order"] = σ_rvs_proj_order
 outputs["ccf_order_resid"] = ccf_order_resid
 outputs["ccf_order_resid_rv_from_order"] = ccf_order_resid_rv_from_order
=#
#=
obs_id = idx_to_fit[6]
 #maxabsz = maximum(abs.(extrema(ccf_order_resid[:,:,obs_id])))
 #plt = heatmap(v_grid_order_ccfs, 1:size(order_ccfs,2), ccf_order_resid[:,:,obs_id]',clims=(-maxabsz,maxabsz),seriescolor=:balance)
 maxabsz = maximum(abs.(extrema(reshape(sum(ccf_order_resid,dims=3),size(ccf_order_resid)[1:2]))))
 plt = heatmap(v_grid_order_ccfs, 1:size(order_ccfs,2), reshape(sum(ccf_order_resid,dims=3),size(ccf_order_resid)[1:2])',clims=(-maxabsz,maxabsz),seriescolor=:balance)
 xlabel!(plt,"v (m/s)")
 ylabel!(plt,"Order (ignore #)")
 title!(plt,"Sun " * date_str  )
 #xlims!(plt,ccf_mid_velocity-RvSpectMLBase.max_bc,ccf_mid_velocity+RvSpectMLBase.max_bc)
 xlims!(plt,ccf_mid_velocity-2.0*line_width_50,ccf_mid_velocity+2.0*line_width_50)
 display(plt)
 #savefig(joinpath(output_dir,"order_ccf_resid_heatmaps_sun_" * date_str * "_espresso_msf=" * string(round(msf,digits=1)) * ".png"))
=#
#=
alg_fit_rv2 =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
order_rvs_g = zeros(length(all_spectra), size(order_ccfs,2))
order_rv_std_g = zeros(length(all_spectra), size(order_ccfs,2))
order_rvs_t = zeros(length(all_spectra), size(order_ccfs,2))
order_rv_std_t = zeros(length(all_spectra), size(order_ccfs,2))
for i in 1:size(order_ccfs,2)
   println("# Index for order: ", i)
   if sum(view(order_ccfs,:,i,1)) > 0
      rvs_order_ccf = RvSpectML.calc_rvs_from_ccf_total(view(order_ccfs,:,i,:), view(order_ccf_vars,:,i,:), pipeline_plan, v_grid=v_grid_order_ccfs, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv2)
      alg_fit_rv3 =  EchelleCCFs.MeasureRvFromCCFTemplate(v_grid=v_grid_order_ccfs, frac_of_width_to_fit=fwtf, template = vec(mean(view(order_ccfs,:,i,:),dims=2)) )
      order_rvs_g[:,i] .= read_cache(pipeline_plan, :rvs_ccf_total )
      order_rv_std_g[:,i] .= read_cache(pipeline_plan, :σ_rvs_ccf_total)
      rvs_order_ccf = RvSpectML.calc_rvs_from_ccf_total(view(order_ccfs,:,i,:), view(order_ccf_vars,:,i,:), pipeline_plan, v_grid=v_grid_order_ccfs, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv3)
      order_rvs_t[:,i] .= read_cache(pipeline_plan, :rvs_ccf_total )
      order_rv_std_t[:,i] .= read_cache(pipeline_plan, :σ_rvs_ccf_total)
   end
end
outputs["order_rvs_g"] = order_rvs_g
outputs["order_rvs_t"] = order_rvs_t
=#

#=
outputs[""] =
outputs[""] =
outputs[""] =
outputs[""] =

=#
#=
order_snr = mapreduce(obsid->map(ord->RvSpectMLBase.calc_snr(all_spectra[obsid],787:6214,ord)  ,1:size(first(all_spectra).flux,2)), hcat, 1:length(all_spectra))
outputs["order_snr"] = order_snr

save(outputs_filename, outputs)
=#
#=

EchelleCCFs.write_each_ccf_fits(map(s->s.metadata,all_spectra), #= 100.0 .* EXPRESS USED cm/s =# v_grid, ccfs_espresso, ccf_vars_espresso,
                  fits_hdr=EchelleCCFs.make_ccf_fits_header(all_spectra[1].metadata),
                  line_list_filename=line_list_filename,# orders=161 .-good_orders,ccf_orders=order_ccfs,ccf_orders_var=order_ccf_vars,
                  total_vel=100.0 .* rvs_ccf_espresso, sigma_total_vel= 100.0 .* σ_rvs_espresso,
                  order_vels = 100.0 .* order_rvs_g, sigma_order_vels = 100.0 .* order_rv_std_g )
=#
