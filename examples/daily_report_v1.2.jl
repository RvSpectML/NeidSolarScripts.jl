verbose = true
 using Dates
 if verbose println(now()) end
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 #using RvSpectMLBase
 #using EchelleInstruments#, EchelleInstruments.NEID
 using RvSpectML
 if verbose println("# Loading NeidSolarScripts")    end
 #using SunAsAStar
 using NeidSolarScripts
 println("# Loading other packages 1/2")
 using ArgParse
 #using Markdown

println("# Parsing arguments...")
   #lsf_width = 3.0e3
   function parse_commandline()
     s = ArgParseSettings( description = "Make dailiy report from daily rvs.")
     #import_settings!(s, s_files_only, args_only=false)
     add_arg_group!(s, "Files", :argg_files)
     @add_arg_table! s begin
         "input"
             help = "Filename for input with daily rvs (csv)"
             arg_type = String
             #default = "daily_rvs.csv"
         "output"
             help = "Filename for output dailiy RVs (toml)"
             arg_type = String
             default = "daily_summary.toml"
         "--overwrite"
            help = "Specify it's ok to overwrite the output file."
            #default = true
            action = :store_true
      end
      add_arg_group!(s, "Filter rvs for ", :argg_filter_param)
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
         "--max_airmass"
            help = "Maximum airmass."
            arg_type = Float64
            default = 2.0
         "--min_expmeter"
            help = "Minimum exposure meter flux relative to model flux."
            arg_type = Float64
            default = 6e4
#=
         "--min_expmeter_factor"
            help = "Minimum exposure meter flux relative to model flux."
            arg_type = Float64
            default = 0.9
         "--max_expmeter_factor"
            help = "Maximum exposure meter flux relative to model flux."
            arg_type = Float64
            default = 1.1
=#
         "--max_expmeter_frac_rms"
            help = "Maximum fractional RMS of exposure meter flux."
            arg_type = Float64
            default = 0.01
#=
         "--min_pyrhelio_factor"
            help = "Minimum pyrheliometer flux relative to model flux."
            arg_type = Float64
            default = 0.9
         "--max_pyrhelio_factor"
            help = "Maximum pyrheliometer flux relative to model flux."
            arg_type = Float64
            default = 1.1
=#
         "--max_pyrhelio_frac_rms"
            help = "Maximum fractional RMS of pyrheliometer meter flux."
            arg_type = Float64
            default = 0.0035
         #=
         "--min_drp_snr"
            help = "Minimum extracted SNR reported by NEID DRP."
            arg_type = Float64
            #default = 0.0
         =#
            #=
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
           =#
         "--max_num_spectra"
            help = "Maximum number of spectra to process."
            arg_type = Int
            default = 300   # TODO: Increase  after finish testing
      end

     return parse_args(s)
 end
 args = parse_commandline()

if verbose   println("# Loading other packages 2/2")    end
 using CSV, DataFrames, Query, Dates
 using FileIO, MD5 #, JLD2, SHA
 using StatsBase, Statistics, NaNMath
 using TOML, OrderedCollections

 # Filename arguments
 @assert isfile(args["input"]) || islink(args["input"])
 @assert match(r"\.csv$",args["input"]) != nothing
 daily_rvs_filename = args["input"]

 args["overwrite"] = true   # TODO: Comment after done testing
 if isfile(args["output"]) && !args["overwrite"]
    error("# Output file " * args["output"] * " already exists (size = " * string(filesize(args["output"])) * " ).")
 end
 @assert !isfile(args["output"]) || args["overwrite"]
 @assert match(r"\.toml$",args["output"]) != nothing
 daily_report_filename = args["output"]
 touch(daily_report_filename)

 file_hashes = Dict{String,String}()
 start_processing_time = now()

println("# Reading daily RVs from ", daily_rvs_filename)
df_in = CSV.read(daily_rvs_filename,DataFrame)

println("# Filtering for usable observations...")

df_use = df_in |>
 #@filter( args["target"] == nothing || _.target == args["target"] ) |>
 @filter( args["datestr"] == nothing || occursin(args["datestr"],_.Filename) ) |>
 #@filter( _.driftfun == "dailymodel0" ) |>
 @filter( args["max_airmass"] == nothing || _.airmass <= args["max_airmass"] ) |>
 @filter( args["max_solar_hour_angle"] == nothing || abs(_.solar_hour_angle) <= args["max_solar_hour_angle"] ) |>
 #@filter( args["start_time"] == nothing || Time(julian2datetime(_.bjd)) >= start_time ) |>
 #@filter( args["stop_time"] == nothing || Time(julian2datetime(_.bjd)) <= stop_time ) |> # TODO for other instruments may need to deal wtih cross end of 24 UTC
 @filter( args["min_expmeter"] == nothing || _.expmeter_mean >= args["min_expmeter"] ) |> 
 DataFrame

if hasproperty(df_use,:pyrflux_mean) && hasproperty(df_use,:pyrflux_rms)
   df_use = df_use |> 
      @filter( args["max_pyrhelio_frac_rms"] == nothing || _.pyrflux_rms <= args["max_pyrhelio_frac_rms"] * _.pyrflux_mean ) |>
      DataFrame
else
   println("# Warning:  mean_pyroflux and rms_pyroflux not found.  Reverting to expmeter_mean and expmeter_rms.")
   df_use = df_use |> 
      @filter( args["max_expmeter_frac_rms"] == nothing || _.expmeter_rms <= args["max_expmeter_frac_rms"] * _.expmeter_mean ) |>
      DataFrame
end

df_use = df_use |> 
   #@take(args["max_num_spectra"] ) |> @orderby(_.bjd) |>
   @take(args["max_num_spectra"] ) |> @orderby(_.jd_drp) |>
   DataFrame

(times_binned, rvs_binned) = bin_times_and_rvs_max_Δt(times=df_use.jd_drp, rvs=df_use.rv_drp, Δt_threshold=6/(60*24))
#rms_rvs_binned = std(rvs_binned,corrected=false)

println("# Found ", size(df_use,1), " files of ",  size(df_use,1), " to use for RVs.")
daily_out = OrderedDict{String,Any}()
#@assert size(df_use,1) >= 1
if size(df_use,1) >= 1
   println("# Creating daily report.")
   daily_out["title"] = "NEID Solar Observations Daily Report"

   obs_date = OrderedDict{String,Any}()
      mean_bjd = mean(df_use.jd_drp)
      obs_date["string"] = isnothing(args["datestr"]) ? Date(julian2datetime(mean_bjd)) : args["datestr"]
      obs_date["mean_bjd"] = mean_bjd
      daily_out["obs_date"] = obs_date
   
   num_rvs = OrderedDict{String,Any}()
     num_rvs["good"] = size(df_use,1)
     num_rvs["usable"] = size(df_in,1)
     daily_out["num_rvs"] = num_rvs

   rv_drp = OrderedDict{String,Any}()
      rv_drp["mean_rv"] = mean(df_use.rv_drp)
      rv_drp["median_rv"] = median(df_use.rv_drp)
      rv_drp["median_σ_rv"] = median(df_use.σrv_drp)
      rv_drp["rms_rvs"] = sqrt(var(df_use.rv_drp,corrected=false))
      rv_drp["rms_binned_rvs"] = sqrt(var(rvs_binned,corrected=false))
      rv_drp["winsor_mean_rv"] = mean(winsor(df_use.rv_drp,prop=0.025))
      rv_drp["winsor_rms_rv"] = sqrt(trimvar(df_use.rv_drp,prop=0.025))
   rv = OrderedDict{String,Any}("drp"=>rv_drp)
   daily_out["rv"] = rv
   
   report = OrderedDict{String,Any}("hostname"=>gethostname(), "date_run"=>now(), "input_file"=>daily_rvs_filename, 
                                "input_md5" => bytes2hex(open(md5,daily_rvs_filename)) )
      daily_out["report"] = report

   println("# Writing daily report.")
   open(daily_report_filename, "w") do io
        TOML.print(io, daily_out)
   end
else
   touch(daily_report_filename)
end


#=
   report_str = """
# NEID Solar Observations Daily Report $(args["datestr"])
- Usable files: $(size(manifest_use,1)) of $(size(manifest,1)) solar observations
- Daily mean JD: $daily_mean_bjd   $(julian2datetime(daily_mean_bjd))
- Daily mean RV: $daily_mean_rv
- Daily median RV: $daily_median_rv
- Median of σ_RVs: $daily_median_σ_rv
- Daily RMS of RVs: $daily_rms_rvs

---
This report generated on $(gethostname()) at $(now()) based on CCFs calculated at $(input_data["start_processing_time"]).
Input file: $(input_data["daily_ccf_filename"])
Input md5sum: $input_md5
"""
else
   report_str = """
# NEID Solar Observations Daily Report $(args["datestr"])
- Usable files: $(size(manifest_use,1)) of $(size(manifest,1)) solar observations

---
This report generated on $(gethostname()) at $(now()) based on CCFs calculated at $(input_data["start_processing_time"]).
Input file: $(input_data["daily_ccf_filename"])
Input md5sum: $input_md5
"""
end

write(daily_report_filename, report_str)
=#

