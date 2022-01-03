verbose = true
 using Dates
 if verbose println(now()) end
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 #using RvSpectMLBase
 #using EchelleInstruments#, EchelleInstruments.NEID
 using RvSpectML
 if verbose println("# Loading NeidSolarScripts")    end
 using NeidSolarScripts
 println("# Loading other packages 1/2")
 using ArgParse
 #using Markdown
 using PkgVersion
 # Just for version info (sigh)
 import EchelleInstruments, EchelleCCFs, RvSpectMLBase, RvSpectML, SunAsAStar

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
#=
         "--template_file"
            help = "Filename with CCF template stored as 'order_ccfs' (jld2) [TODO]"
            arg_type = String
=#
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
warning_str = ""

println("# Filtering for usable observations...")

df_use = df_in |>
 #@filter( args["target"] == nothing || _.target == args["target"] ) |>
 @filter( args["datestr"] == nothing || occursin(args["datestr"],_.Filename) ) |>
 #@filter( _.driftfun == "dailymodel0" ) |>
 @filter( args["max_airmass"] == nothing || _.airmass <= args["max_airmass"] ) |>
 @filter( args["max_solar_hour_angle"] == nothing || abs(_.hour_angle) <= args["max_solar_hour_angle"] ) |>
 #@filter( args["start_time"] == nothing || Time(julian2datetime(_.bjd)) >= start_time ) |>
 #@filter( args["stop_time"] == nothing || Time(julian2datetime(_.bjd)) <= stop_time ) |> # TODO for other instruments may need to deal wtih cross end of 24 UTC
 @filter( args["min_expmeter"] == nothing || _.expmeter_mean >= args["min_expmeter"] ) |> 
 DataFrame

if hasproperty(df_use,:pyrflux_mean) && hasproperty(df_use,:pyrflux_rms)
   df_use = df_use |> 
      @filter( args["max_pyrhelio_frac_rms"] == nothing || _.pyrflux_rms <= args["max_pyrhelio_frac_rms"] * _.pyrflux_mean ) |>
      DataFrame
else
   add_warning_str = "" * "# Warning:  mean_pyroflux and rms_pyroflux not found.  Reverting to expmeter_mean and expmeter_rms."
   warning_str = "" * add_warning_str
   println(add_warning_str)
   df_use = df_use |> 
      @filter( args["max_expmeter_frac_rms"] == nothing || _.expmeter_rms <= args["max_expmeter_frac_rms"] * _.expmeter_mean ) |>
      DataFrame
end

df_use = df_use |> 
   #@take(args["max_num_spectra"] ) |> @orderby(_.bjd) |>
   @take(args["max_num_spectra"] ) |> @orderby(_.jd_drp) |>
   DataFrame


println("# Found ", size(df_use,1), " files of ",  size(df_use,1), " to use for RVs.")
daily_out = OrderedDict{String,Any}()
#@assert size(df_use,1) >= 1
if size(df_use,1) >= 1
   println("# Creating daily report.")
   daily_out["title"] = "NEID Solar Observations Daily Report"
   if length(warning_str) > 0
      daily_out["warning"] = warning_str
   end

   obs_date = OrderedDict{String,Any}()
      mean_bjd = mean(df_use.jd_drp)
      obs_date["string"] = isnothing(args["datestr"]) ? Date(julian2datetime(mean_bjd)) : args["datestr"]
      obs_date["mean_bjd"] = mean_bjd
      daily_out["obs_date"] = obs_date
   
   num_rvs = OrderedDict{String,Any}()
     num_rvs["good"] = size(df_use,1)
     num_rvs["usable"] = size(df_in,1)
     daily_out["num_rvs"] = num_rvs

   rv = OrderedDict{String,Any}()
   rv_drp = OrderedDict{String,Any}()
      drp_vel_to_mps = 1000
      rv_drp["mean_rv"] = mean(df_use.rv_drp)  * drp_vel_to_mps
      rv_drp["median_rv"] = median(df_use.rv_drp) * drp_vel_to_mps
      rv_drp["median_σ_rv"] = median(df_use.σrv_drp) * drp_vel_to_mps
      rv_drp["rms_rvs"] = sqrt(var(df_use.rv_drp,corrected=false)) * drp_vel_to_mps
      (times_binned, rvs_binned) = bin_times_and_rvs_max_Δt(times=df_use.jd_drp, rvs=df_use.rv_drp, Δt_threshold=6/(60*24))
      rv_drp["rms_binned_rvs"] = sqrt(var(rvs_binned,corrected=false)) * drp_vel_to_mps
      
      trimed_rvs = trim(df_use.rv_drp,prop=0.025) .* drp_vel_to_mps
      if length(collect(trimed_rvs)) >= 2
      rv_drp["winsor_mean_rv"] = mean(trimed_rvs)
      rv_drp["winsor_rms_rv"] = sqrt(mean((trimed_rvs.- rv_drp["winsor_mean_rv"]).^2))
      end
      trimed_binned_rvs = trim(rvs_binned,prop=0.025) .* drp_vel_to_mps
      if length(collect(trimed_binned_rvs)) >= 2
      mean_binned_rvs = mean(trimed_binned_rvs)
      rv_drp["winsor_rms_binned_rvs"] = sqrt(mean((trimed_binned_rvs.- mean_binned_rvs).^2))
      end
   rv["drp"] = rv_drp

   if hasproperty(df_use,:rv_template)
   rv_template= OrderedDict{String,Any}()
      rv_template["mean_rv"] = mean(df_use.rv_template)
      rv_template["median_rv"] = median(df_use.rv_template)
      rv_template["median_σ_rv"] = median(df_use.σrv_template)
      rv_template["rms_rvs"] = sqrt(var(df_use.rv_template,corrected=false))
      (times_binned, rvs_binned) = bin_times_and_rvs_max_Δt(times=df_use.jd_drp, rvs=df_use.rv_template, Δt_threshold=6/(60*24))
      rv_template["rms_binned_rvs"] = sqrt(var(rvs_binned,corrected=false))
      
      trimed_rvs = trim(df_use.rv_template,prop=0.025)
      if length(collect(trimed_rvs)) >= 2
      rv_template["winsor_mean_rv"] = mean(trimed_rvs)
      rv_template["winsor_rms_rv"] = sqrt(mean((trimed_rvs.- rv_template["winsor_mean_rv"]).^2))
      end
      trimed_binned_rvs = trim(rvs_binned,prop=0.025)
      if length(collect(trimed_binned_rvs)) >= 2
      mean_binned_rvs = mean(trimed_binned_rvs)
      rv_template["winsor_rms_binned_rvs"] = sqrt(mean((trimed_binned_rvs.- mean_binned_rvs).^2))
      end
      rv["template"] = rv_template
   end

   if hasproperty(df_use,:rv_gaussian)
   rv_gauss= OrderedDict{String,Any}()
      rv_gauss["mean_rv"] = mean(df_use.rv_gaussian)
      rv_gauss["median_rv"] = median(df_use.rv_gaussian)
      rv_gauss["median_σ_rv"] = median(df_use.σrv_gaussian)
      rv_gauss["rms_rvs"] = sqrt(var(df_use.rv_gaussian,corrected=false))
      (times_binned, rvs_binned) = bin_times_and_rvs_max_Δt(times=df_use.jd_drp, rvs=df_use.rv_gaussian, Δt_threshold=6/(60*24))
      rv_gauss["rms_binned_rvs"] = sqrt(var(rvs_binned,corrected=false))
      
      trimed_rvs = trim(df_use.rv_gaussian,prop=0.025)
      if length(collect(trimed_rvs)) >= 2
      rv_gauss["winsor_mean_rv"] = mean(trimed_rvs)
      rv_gauss["winsor_rms_rv"] = sqrt(mean((trimed_rvs.- rv_gauss["winsor_mean_rv"]).^2))
      end
      trimed_binned_rvs = trim(rvs_binned,prop=0.025)
      if length(collect(trimed_binned_rvs)) >= 2
      mean_binned_rvs = mean(trimed_binned_rvs)
      rv_gauss["winsor_rms_binned_rvs"] = sqrt(mean((trimed_binned_rvs.- mean_binned_rvs).^2))
      end
      rv["gaussian"] = rv_gauss
   end

   #rv = OrderedDict{String,Any}("drp"=>rv_drp, "ccf_template"=>rv_template, "ccf_gaussian"=>rv_gauss )
   daily_out["rv"] = rv
   
   if hasproperty(df_use, :CaIIHK) && hasproperty(df_use, :σ_CaIIHK)
      indicator = OrderedDict{String,Any}()
      indicator["CaIIHK"] = mean(df_use.CaIIHK)
      indicator["σ_CaIIHK"] = median(df_use.σ_CaIIHK)
      indicator["Ha06_1"] = mean(df_use.Ha06_1)
      indicator["σ_Ha06_1"] = median(df_use.σ_Ha06_1)
      indicator["Ha06_2"] = mean(df_use.Ha06_2)
      indicator["σ_Ha06_2"] = median(df_use.σ_Ha06_2)
      indicator["Ha16_1"] = mean(df_use.Ha16_1)
      indicator["σ_Ha16_1"] = median(df_use.σ_Ha16_1)
      indicator["Ha16_2"] = mean(df_use.Ha16_2)
      indicator["σ_Ha16_2"] = median(df_use.σ_Ha16_2)
      #indicator["NaI"] = df_use.NaI
      #indicator["σ_NaI"] = df_use.σ_NaI
      daily_out["indicators"] = indicator 
   end

   quality = OrderedDict{String,Any}()
      quality["extsnr"] = median(df_use.drpextsnr)
      if hasproperty(df_use, :expmeter_mean) && hasproperty(df_use, :pyrflux_mean) 
      quality["median_expmeter_over_pyrheliometer"] = median(df_use.expmeter_mean ./ df_use.pyrflux_mean)
      quality["min_expmeter_over_pyrheliometer"] = minimum(df_use.expmeter_mean ./ df_use.pyrflux_mean)
      quality["max_expmeter_over_pyrheliometer"] = maximum(df_use.expmeter_mean ./ df_use.pyrflux_mean)
      end
   daily_out["quality"] = quality
   environment = OrderedDict{String,Any}()
      if hasproperty(df_use, :envohum)   environment["humidity"] = median(df_use.envohum) end
      if hasproperty(df_use, :envotmp)   environment["temperature"] = median(df_use.envotmp) end
      if hasproperty(df_use, :envwinds)  environment["wind_speed"] = median(df_use.envwinds) end
      if hasproperty(df_use, :envwindd)  environment["wind_direction"] = median(df_use.envwindd) end
   daily_out["environment"] = environment

   report = OrderedDict{String,Any}()
      report["hostname"]=gethostname()
      report["date_run"]=now() 
      report["input_file"]=daily_rvs_filename
      report["input_ctime"] = ctime(daily_rvs_filename)
      report["input_md5"] = bytes2hex(open(md5,daily_rvs_filename)) 
      report_versions = OrderedDict{String,Any}()
      report_versions["NeidSolarScripts"] = string(PkgVersion.Version(NeidSolarScripts))
      report_versions["EchelleCCFs"] = string(PkgVersion.Version(EchelleCCFs))
      report_versions["EchelleInstruments"] = string(PkgVersion.Version(EchelleInstruments))
      report_versions["RvSpectML"] = string(PkgVersion.Version(RvSpectML))
      report_versions["RvSpectMLBase"] = string(PkgVersion.Version(RvSpectMLBase))
      report_versions["SunAsAStar"] = string(PkgVersion.Version(SunAsAStar))
      report["versions"] = report_versions
   daily_out["report"] = report

   println("# Writing daily report.")
   open(daily_report_filename, "w") do io
        TOML.print(io, daily_out) end
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

