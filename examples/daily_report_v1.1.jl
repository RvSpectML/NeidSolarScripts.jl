verbose = true
 using Dates
 if verbose println(now()) end
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectMLBase
 using EchelleInstruments#, EchelleInstruments.NEID
 using RvSpectML
 if verbose println("# Loading NeidSolarScripts")    end
 using SunAsAStar
 using NeidSolarScripts
 println("# Loading other packages 1/2")
 using ArgParse
 using Markdown

println("# Parsing arguments...")
   #lsf_width = 3.0e3
   function parse_commandline()
     s = ArgParseSettings( description = "Calculate RVs from JLD2 file with NEID's daily order CCFs.")
     #import_settings!(s, s_files_only, args_only=false)
     add_arg_group!(s, "Files", :argg_files)
     @add_arg_table! s begin
         "input"
             help = "Filename for input with daily order CCFs (jld2)"
             arg_type = String
             default = "daily_ccfs.jld2"
         "output"
             help = "Filename for output dailiy RVs (csv)"
             arg_type = String
             default = "daily_rvs.csv"
         "report"
             help = "Filename for summary repork (md)"
             arg_type = String
             default = "daily_summary.md"
         "--overwrite"
            help = "Specify it's ok to overwrite the output file."
            #default = true
            action = :store_true
      end
      #=
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
      end
      =#
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
           default = 6.0
         "--max_airmass"
            help = "Maximum airmass."
            arg_type = Float64
            default = 10.0
         "--min_expmeter_factor"
            help = "Minimum exposure meter flux relative to model flux."
            arg_type = Float64
            default = 0.9
         "--max_expmeter_factor"
            help = "Maximum exposure meter flux relative to model flux."
            arg_type = Float64
            default = 1.1
         "--max_expmeter_frac_rms"
            help = "Maximum fractional RMS of exposure meter flux."
            arg_type = Float64
            default = 0.01
         "--min_pyrohelio_factor"
            help = "Minimum pyroheliometer flux relative to model flux."
            arg_type = Float64
            default = 0.9
         "--max_pyrohelio_factor"
            help = "Maximum pyroheliometer flux relative to model flux."
            arg_type = Float64
            default = 1.1
         "--max_pyrohelio_frac_rms"
            help = "Maximum fractional RMS of pyroheliometer meter flux."
            arg_type = Float64
            default = 0.01
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
 using JLD2, FileIO, MD5 #, SHA
 using StatsBase, Statistics, NaNMath

 # Filename arguments
 @assert isfile(args["input"]) || islink(args["input"])
 @assert match(r"\.jld2$",args["input"]) != nothing
 daily_ccf_filename = args["input"]

 args["overwrite"] = true   # TODO: Comment after done testing
 if isfile(args["output"]) && !args["overwrite"]
    error("# Output file " * args["output"] * " already exists (size = " * string(filesize(args["output"])) * " ).")
 end
 @assert !isfile(args["output"]) || args["overwrite"]
 @assert match(r"\.csv$",args["output"]) != nothing
 daily_rvs_filename = args["output"]
 touch(daily_rvs_filename)

 if isfile(args["report"]) && !args["overwrite"]
    error("# Report file " * args["report"] * " already exists (size = " * string(filesize(args["report"])) * " ).")
 end
 @assert !isfile(args["report"]) || args["overwrite"] == true
 @assert match(r"\.md$",args["report"]) != nothing
 daily_report_filename = args["report"]
 touch(daily_report_filename)

 file_hashes = Dict{String,String}()
 start_processing_time = now()

println("# Reading daily CCFs from ", daily_ccf_filename)
input_data = load(daily_ccf_filename)

println("# Filtering for usable observations...")
manifest = input_data["manifest"]

manifest_use = manifest |>
 @filter( args["target"] == nothing || _.target == args["target"] ) |>
 @filter( args["datestr"] == nothing || occursin(args["datestr"],_.target) ) |>
 @filter( _.driftfun == "dailymodel0" ) |>
 @filter( args["max_airmass"] == nothing || _.airmass <= args["max_airmass"] ) |>
 @filter( args["max_solar_hour_angle"] == nothing || abs(_.solar_hour_angle) <= args["max_solar_hour_angle"] ) |>
 #@filter( args["start_time"] == nothing || Time(julian2datetime(_.bjd)) >= start_time ) |>
 #@filter( args["stop_time"] == nothing || Time(julian2datetime(_.bjd)) <= stop_time ) |> # TODO for other instruments may need to deal wtih cross end of 24 UTC
 # TODO: Replace to use expmeter or pyrheliometer data
 # @filter( args["min_snr_factor"] == nothing || sum(_.order_snrs) >= args["min_snr_factor"] * max_snr ) |>
 @filter( args["max_expmeter_frac_rms"] == nothing || _.expmeter_rms <= args["max_expmeter_frac_rms"] * _.expmeter_mean ) |>
 @take(args["max_num_spectra"] ) |> @orderby(_.bjd) |>
 DataFrame
println("# Found ", size(manifest_use,1), " files of ",  size(manifest,1), " to use for RVs.")
@assert size(manifest_use,1) >= 1

df_out = select(manifest_use,[:drp_ccfjdmod=>:jd_drp,:drp_ccfrvmod => :rv_drp,:drp_dvrmsmod => :σrv_drp ], :Δv_diff_ext, :Δfwhm², :solar_hour_angle, :airmass, :sol_dist, :expmeter_mean, :expmeter_rms, [:mean_pyroflux=>:pyroflux_mean, :rms_pyroflux=>:pyroflux_rms], :exptime, :mean_Δt, :Filename)

daily_mean_bjd = mean(df_out.jd_drp)
daily_mean_rv = mean(df_out.rv_drp)
daily_median_rv = median(df_out.rv_drp)
daily_median_σ_rv = median(df_out.σrv_drp)
daily_rms_rvs = sqrt(var(df_out.rv_drp,corrected=false))

println("# Writing outputs")
CSV.write(daily_rvs_filename,df_out)

input_md5 = bytes2hex(open(md5,daily_ccf_filename))

report_str = """
# NEID Solar Observations Daily Report
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

write(daily_report_filename, report_str)
