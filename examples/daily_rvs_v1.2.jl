verbose = true
 using Dates
 if verbose println(now()) end
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectMLBase
 using EchelleCCFs
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
      add_arg_group!(s, "Algorithm more measuring RVs ", :argg_alg_param)
      @add_arg_table! s begin
         "--rv_alg"
            help = "Algorithm for measuring RVs."
            arg_type = String
            default = "template"
            #default = "gaussian"
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
 @assert in(args["rv_alg"], ["template","gaussian","quadratic"] )

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

#=
 if isfile(args["report"]) && !args["overwrite"]
    error("# Report file " * args["report"] * " already exists (size = " * string(filesize(args["report"])) * " ).")
 end
 @assert !isfile(args["report"]) || args["overwrite"] == true
 @assert match(r"\.md$",args["report"]) != nothing
 daily_report_filename = args["report"]
 touch(daily_report_filename)
=#

 file_hashes = Dict{String,String}()
 start_processing_time = now()

if filesize(daily_ccf_filename) > 0
   println("# Reading daily CCFs from ", daily_ccf_filename)
   input_data = load(daily_ccf_filename)
else
   println("# Empy daily CCF file.  Creating empty ", daily_rvs_filename)
   touch(daily_rvs_filename)
   exit(0)
end

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
 @filter( args["min_expmeter"] == nothing || _.expmeter_mean >= args["min_expmeter"] ) |> 
 DataFrame

if hasproperty(manifest_use,:mean_pyroflux) && hasproperty(manifest_use,:rms_pyroflux)
   manifest_use = manifest_use |> 
      @filter( args["max_pyrhelio_frac_rms"] == nothing || _.rms_pyroflux <= args["max_pyrhelio_frac_rms"] * _.mean_pyroflux ) |>
      DataFrame
else
   println("# Warning:  mean_pyroflux and rms_pyroflux not found.  Reverting to expmeter_mean and expmeter_rms.")
   manifest_use = manifest_use |> 
      @filter( args["max_expmeter_frac_rms"] == nothing || _.expmeter_rms <= args["max_expmeter_frac_rms"] * _.expmeter_mean ) |>
      DataFrame
end

manifest_use = manifest_use |> 
   @take(args["max_num_spectra"] ) |> @orderby(_.bjd) |>
   DataFrame

#@assert size(manifest_use,1) >= 1
if !(size(manifest_use,1) >= 1)
   println("# No usable observations.  Creating empty ", daily_rvs_filename)
   touch(daily_rvs_filename)
   exit(0)
end

println("# Found ", size(manifest_use,1), " files of ",  size(manifest,1), " to use for RVs.")
   #df_out = manifest_use |> @rename(:drp_ccfjdmod=>:jd_drp,:drp_ccfrvmod => :rv_drp,:drp_dvrmsmod => :σrv_drp, 
   df_out = manifest |> @rename(:drp_ccfjdmod=>:jd_drp,:drp_ccfrvmod => :rv_drp,:drp_dvrmsmod => :σrv_drp, 
                                    :mean_pyroflux=>:pyrflux_mean, :rms_pyroflux=>:pyrflux_rms ) |> DataFrame 
   obs_to_use = map(i->searchsortedfirst(manifest.bjd,manifest_use.bjd[i]),1:size(manifest_use,1) )

(num_vels, num_orders_to_use, num_obs) = size(input_data["order_ccfs"])
 
println("# Measuring RVs from RvSpectML's CCFs.")
v_grid = collect(input_data["v_grid"])
@assert num_vels == length(v_grid)
order_ccfs = input_data["order_ccfs"]
non_nan_mask = map(i->!any(isnan.(view(order_ccfs,:,:,i))), obs_to_use)
obs_to_use = obs_to_use[non_nan_mask]
   mask = falses(1:num_obs)
   mask[obs_to_use] .= true
   df_out[!,"mask"] = mask

order_ccf_vars = input_data["order_ccf_vars"]
ccf_sum = reshape(sum(order_ccfs,dims=2),num_vels,num_obs)
ccf_var_sum = reshape(sum(input_data["order_ccf_vars"],dims=2),num_vels,num_obs)
#rvs = measure_rv_from_ccf(v_grid,view(ccf_sum,:,1),view(ccf_var_sum,:,1))
if args["rv_alg"] == "template"
   ccf_template = EchelleCCFs.calc_ccf_template(view(ccf_sum,:,obs_to_use), view(ccf_var_sum,:,obs_to_use) )
   alg_template = MeasureRvFromCCFTemplate(v_grid=collect(v_grid),template=ccf_template)
   rvs_template = DataFrame(map(i->measure_rv_from_ccf(v_grid,view(ccf_sum,:,i),view(ccf_var_sum,:,i),alg=alg_template),1:num_obs))
   df_out[!,Symbol("rv_template")] = rvs_template.rv
   df_out[!,Symbol("σrv_template")] = rvs_template.σ_rv
elseif args["rv_alg"] == "gaussian"
   alg_gauss = MeasureRvFromCCFGaussian(frac_of_width_to_fit=0.25)
   rvs_gauss = DataFrame(map(i->measure_rv_from_ccf(v_grid,view(ccf_sum,:,i),view(ccf_var_sum,:,i),alg=alg_gauss),1:num_obs))
   df_out[!,Symbol("rv_gaussian")] = rvs_gauss.rv
   df_out[!,Symbol("σrv_gaussian")] = rvs_gauss.σ_rv
elseif args["rv_alg"] == "quadratic"
   alg_quad = MeasureRvFromCCFQuadratic(frac_of_width_to_fit=0.25)
   rvs_quad = DataFrame(map(i->measure_rv_from_ccf(v_grid,view(ccf_sum,:,i),view(ccf_var_sum,:,i),alg=alg_quad),1:num_obs))
   df_out[!,Symbol("rv_quadratic")] = rvs_quad.rv
   df_out[!,Symbol("σrv_quadratic")] = rvs_quad.σ_rv
end

orders_physical = collect(174 .-input_data["orders_to_use"])
for j in 1:num_orders_to_use
   if iszero(maximum(view(order_ccfs,:,j,:))) || any(isnan.(view(order_ccfs,:,j,:)))   continue    end
   #println(" j = ", j, "  order_fits_idx = ", input_data["orders_to_use"][j], "  order_physical = ", orders_physical[j])
   if args["rv_alg"] == "template"
   # Templated-based
   local ccf_template = EchelleCCFs.calc_ccf_template(view(order_ccfs,:,j,obs_to_use), view(order_ccf_vars,:,j,obs_to_use) )
   local alg_template = MeasureRvFromCCFTemplate(v_grid=collect(v_grid),template=ccf_template)
   order_rvs_template = DataFrame(map(i->measure_rv_from_ccf(v_grid,view(order_ccfs,:,j,i),view(order_ccf_vars,:,j,i),alg=alg_template),1:num_obs)) 
   df_out[!,Symbol("rv_" * string(orders_physical[j]) * "_template")] = order_rvs_template.rv
   df_out[!,Symbol("σrv_" * string(orders_physical[j]) * "_template")] = order_rvs_template.σ_rv
   elseif args["rv_alg"] == "gaussian"
   # Gaussian fitting
   order_rvs_gauss = DataFrame(map(i->measure_rv_from_ccf(v_grid,view(order_ccfs,:,j,i),view(order_ccf_vars,:,j,i),alg=alg_gauss),1:num_obs))
   df_out[!,Symbol("rv_" * string(orders_physical[j]) * "_gauss")] = order_rvs_gauss.rv
   df_out[!,Symbol("σrv_" * string(orders_physical[j]) * "_gauss")] = order_rvs_gauss.σ_rv
   elseif args["rv_alg"] == "quadratic"
   # Quadratic fitting
   order_rvs_quad = DataFrame(map(i->measure_rv_from_ccf(v_grid,view(order_ccfs,:,j,i),view(order_ccf_vars,:,j,i),alg=alg_quad),1:num_obs))
   df_out[!,Symbol("rv_" * string(orders_physical[j]) * "_quad")] = order_rvs_quad.rv
   df_out[!,Symbol("σrv_" * string(orders_physical[j]) * "_quad")] = order_rvs_quad.σ_rv
   end
end

println("# Writing daily RVs.")
 CSV.write(daily_rvs_filename,df_out)

