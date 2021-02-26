if occursin(r"RvSpectMLEcoSystem$", pwd())   cd("NeidSolarScripts")   end
 using Pkg
 Pkg.activate(".")

verbose = true
 make_plots = true
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectMLBase
 using EchelleInstruments, EchelleInstruments.NEID
 using EchelleCCFs
 using RvSpectML
 if verbose   println("# Loading other packages")    end
 using CSV, DataFrames, Query, StatsBase, Statistics, Dates
 using JLD2, FileIO
 using NaNMath

target_subdir = "good_days"   # USER: Replace with directory of your choice
  fits_target_str = "Sun"
  output_dir = "output/"
  outputs = Dict{String,Any}()
  paths_to_search_for_param = [pwd(),"examples",joinpath(pkgdir(RvSpectMLBase),"..","RvSpectML","examples"), "/gpfs/scratch/jpn23/"]
  # NOTE: make_manifest does not update its paths_to_search when default_paths_to_search is defined here, so if you change the line above, you must also include "paths_to_search=default_paths_to_search" in the make_manifest() function call below
  pipeline_plan = PipelinePlan()
  dont_make_plot!(pipeline_plan, :movie)

reset_all_needs!(pipeline_plan)
#if need_to(pipeline_plan,:read_spectra)
if verbose println("# Finding what data files are avaliable.")  end
if isfile("manifest.csv")
    df_files  = CSV.read("manifest.csv", DataFrame)
    @assert size(df_files,1) >= 1
    @assert hasproperty(df_files,:Filename)
    @assert hasproperty(df_files,:target)
    @assert hasproperty(df_files,:bjd)
    @assert hasproperty(df_files,:ssbz)
    @assert hasproperty(df_files,:exptime)
else
    eval(read_data_paths(paths_to_search=paths_to_search_for_param))
    @assert isdefined(Main,:neid_data_path)
    df_files = make_manifest(neid_data_path, target_subdir, NEID )
    CSV.write("manifest.csv", df_files)
end

if verbose println("# Reading in customized parameters from param.jl.")  end
   if !@isdefined(idx_day_to_use)
       idx_day_to_use = 1
   end
   eval(code_to_include_param_jl(paths_to_search=paths_to_search_for_param))
   if match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1] ==  match(r"neidL1_(\d+)[T_](\d+)\.fits$", last(df_files_use.Filename))[1]
      date_str = match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1]
    else
      date_str = string(match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1]) * "-" * string(match(r"neidL1_(\d+)[T_](\d+)\.fits$", last(df_files_use.Filename))[1])
   end
   outputs["df_files_use"] = df_files_use

   outputs_filename = joinpath(output_dir,"solar_" * date_str * ".jld2")
   if isfile(outputs_filename) && false
     times_already_processed = load(outputs_filename, "times")
     files_in_day_to_process = size(df_files_solar_by_day.data[idx_day_to_use],1)
      if files_in_day_to_process == length(times_already_processed)
         println("# Already processed all ", length(times_already_processed), " files for ", date_str)
         exit()
      end
   end

#for idx_day_to_use in 1:size(df_files_solar_by_day,1)
   #df_files_use = df_files_solar_by_day[idx_day_to_use,:data] |> @orderby(_.bjd) |> @take(max_spectra_to_use) |> DataFrame

if verbose println("# Reading in ", size(df_files_use,1), " FITS files for ", date_str, ".")  end
    @time all_spectra = map(row->NEID.read_data(row), eachrow(df_files_use))
    GC.gc()
    dont_need_to!(pipeline_plan,:read_spectra)


line_width_50_default = 7.9e3
 lsf_width = 3.0e3
 max_mask_scale_factor = 4.0
 max_bc = RvSpectMLBase.max_bc
 #max_bc = RvSpectMLBase.max_bc_earth_rotation
 if isfile(joinpath("data","solar_line_list_espresso.csv"))
   line_list_espresso = CSV.read(joinpath("data","solar_line_list_espresso.csv"), DataFrame)
   dont_need_to!(pipeline_plan,:clean_line_list_tellurics)
 else
    #max_orders = min_order(NEID2D()):max_order(NEID2D())
    good_orders = orders_to_use_default(NEID2D())
    #order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=max_orders, recalc=true )
    order_list_timeseries = extract_orders(all_spectra,pipeline_plan,  orders_to_use=good_orders,  remove_bad_chunks=false, recalc=true )

    line_list_filename = "G2.espresso.mas"
    linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",line_list_filename)

    line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity,
       Δv_to_avoid_tellurics = 2*max_bc+5*line_width_50_default+max_mask_scale_factor*default_ccf_mask_v_width(NEID1D()), orders_to_use=good_orders, recalc=true, verbose=true)

    #CSV.write(joinpath("data","solar_line_list_espresso.csv"), line_list_espresso)
 end
 outputs["line_list_espresso"] = line_list_espresso
 line_width_50 = line_width_50_default



good_orders = 56:111 # orders_to_use_default(NEID2D())
 #good_orders = 57:96 # orders_to_use_default(NEID2D())
 #good_orders = 55 .+ vcat(1:14, 17:24, 27:29, 37:38,40:41)
 order_list_timeseries = extract_orders(all_spectra,pipeline_plan,  orders_to_use=good_orders, remove_bad_chunks=false, recalc=true )
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

for obsid in 1:length(order_list_timeseries)
    for ch in 1:length(order_list_timeseries[1].data)
        order_list_timeseries[obsid].data[ch].flux .*= order_weights[ch,obsid] #./mean_order_mean_flux[ch]
        order_list_timeseries[obsid].data[ch].var .*= (order_weights[ch,obsid]).^2 # ./mean_order_mean_flux[ch]).^2
    end
end

#=for obsid in 1:length(order_list_timeseries)
    for ch in 1:length(order_list_timeseries[1].data)
        order_list_timeseries[obsid].data[ch].flux ./= order_weights[ch,obsid] #./mean_order_mean_flux[ch]
        order_list_timeseries[obsid].data[ch].var ./= (order_weights[ch,obsid]).^2 # ./mean_order_mean_flux[ch]).^2
    end
end
=#

using Plots
 plot(map(obsid->NaNMath.sum(order_list_timeseries[obsid].data[10].flux)/(NaNMath.sum(order_list_timeseries[obsid].data[30].flux)), 1:length(all_spectra)),label="Raw")
 #plot!(map(obsid->NaNMath.sum(order_list_timeseries[obsid].data[10].flux)*order_weights[10,obsid]/(NaNMath.sum(order_list_timeseries[obsid].data[30].flux*order_weights[30,obsid])), 1:length(all_spectra)),label="weighted")
 #plot!(map(obsid->NaNMath.sum(order_list_timeseries[obsid].data[10].flux)*order_weights[10,obsid]^2/(NaNMath.sum(order_list_timeseries[obsid].data[30].flux*order_weights[30,obsid]^2)), 1:length(all_spectra)),label="weighted^2")

# plot(order_list_timeseries[obsid].data[ch].λ,order_list_timeseries[obsid].data[ch].flux)


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

#msf = 1.4; fwtf = 1.5
#msf = 6.0; fwtf = 1.5
msf = 0.62; fwtf = 1.5  # one pixel because tophat not using NEID's default v width
#msf = 1.24; fwtf = 1.5   # two pixels
 println("# mask_scale_factor = ", msf, ", frac_width_to_fit =  ", fwtf, ". ")
 ((ccfs_espresso, ccf_vars_espresso), v_grid) = ccf_total(order_list_timeseries, line_list_espresso, pipeline_plan,
   #mask_type=:gaussian,
   mask_scale_factor=msf, range_no_mask_change=5*line_width_50, ccf_mid_velocity=ccf_mid_velocity, v_step=155,
   v_max=max(5*line_width_50,2*max_bc), allow_nans=true, calc_ccf_var=true, recalc=true)
  outputs["v_grid"] = v_grid
  outputs["ccfs_espresso"] = ccfs_espresso
  outputs["ccf_vars_espresso"] = ccf_vars_espresso

fwtf = 1.5; println("# mask_scale_factor = ", msf, ", frac_width_to_fit =  ", fwtf, ". ")
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
 times .-= time_max_alt_sun
 #times = order_list_timeseries.times.-first(order_list_timeseries.times)
 x = [times[idx_to_fit]  ones(length(times[idx_to_fit])) ]
 linfit = (x'*x)\x'*(rvs_ccf_espresso[idx_to_fit])
 x = [times  ones(length(times)) ]
 v_fit = (x*linfit)
 v_t0 = [[0] [1]]*linfit
 (times_binned, rvs_binned) = bin_times_and_rvs_max_Δt(times=times, rvs=rvs_ccf_espresso, Δt_threshold=5/(60*24))
 #idx_to_fit_binned = findfirst(t->t>=times[first(idx_to_fit)], times_binned):findlast(t->t<=times[last(idx_to_fit)], times_binned)
 idx_to_fit_binned = findfirst(t->time_max_alt_sun-t<=Δt_fit,times_binned):findlast(t->t-time_max_alt_sun<=Δt_fit,times_binned)
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
 plt = scatter( (times).*24, rvs_ccf_espresso.-v_fit,  ms=1.5, label=:none)
 scatter!(plt,(times_binned).*24, rvs_binned.-v_fit_binned, ms=4.0, label=:none)
 vline!(plt,[-Δt_fit,Δt_fit].* 24, label=:none, color=:grey)
 vline!(plt,[0].* 24, label=:none, color=:orange)
 xlims!(plt,-3.5,3.5)
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
 maxabsz = maximum(abs.(extrema(view(ccfs_norm,cols_to_fit,:)'.-view(ccf_template,cols_to_fit,:)')))
 plt = heatmap(view(v_grid,cols_to_fit),(times[1:size(ccfs_norm,2)]).*24,view(ccfs_norm,cols_to_fit,:)'.-view(ccf_template,cols_to_fit,:)',clims=(-maxabsz,maxabsz),seriescolor=:balance)
 xlabel!(plt,"v (m/s)")
 ylabel!(plt,"Observation Time (hr)")
 title!(plt,"Sun " * date_str * ":\nCCF - <CCF>, ESPRESSO mask + SNR weights")
 #xlims!(ccf_mid_velocity-RvSpectMLBase.max_bc,ccf_mid_velocity+RvSpectMLBase.max_bc)
 xlims!(plt,ccf_mid_velocity-2*line_width_50,ccf_mid_velocity+2*line_width_50)
 hline!(plt,[-Δt_fit,Δt_fit].*24, label=:none, color=:grey)
 hline!(plt,[0.], label=:none, color=:orange)
 ylims!(plt,-3.5,3.5)
 #vline!([-2.1e4], color=4, label=:none)
 display(plt)
 savefig(joinpath(output_dir,"ccf_heatmaps_sun_" * date_str * "_espresso_msf=" * string(round(msf,digits=1)) * ".png"))

 outputs["mean_ccf"] = ccf_template

gp_post = RvSpectML.TemporalGPInterpolation.construct_gp_posterior(v_grid,ccf_template,sigmasq_obs=ccf_vars_template,use_logx=false,smooth_factor=2)
 ccf_template_smooth = RvSpectML.TemporalGPInterpolation.predict_mean(gp_post(v_grid))
 ccf_template_deriv = RvSpectML.numerical_deriv(v_grid,RvSpectML.TemporalGPInterpolation.predict_mean(gp_post(v_grid)))
 norm = sum(abs2.(ccf_template_deriv[cols_to_fit]))
 rvs_proj = vec(sum((ccf_template_smooth[cols_to_fit].-ccfs_norm[cols_to_fit,:]).*ccf_template_deriv[cols_to_fit],dims=1)./norm)
 σ_rvs_proj = vec(sum((ccf_vars_norm[cols_to_fit]).*abs2.(ccf_template_deriv[cols_to_fit]),dims=1))./ norm
 ccf_resid_minus_rv_proj = (ccfs_norm.-ccf_template).+(rvs_proj.*ccf_template_deriv')'

 outputs["ccf_template_smooth"] = ccf_template_smooth
 outputs["rvs_proj"] = rvs_proj
 outputs["σ_rvs_proj"] = σ_rvs_proj
 outputs["ccf_resid_minus_rv_proj"] = ccf_resid_minus_rv_proj

#plt = heatmap(view(v_grid,cols_to_fit),(times[1:size(ccfs_norm,2)]).*24,ccf_resid_minus_rv_proj[cols_to_fit,:]')
maxabsz = maximum(abs.(extrema(ccf_resid_minus_rv_proj)))
 plt = heatmap(v_grid,(times[1:size(ccfs_norm,2)]).*24,ccf_resid_minus_rv_proj',clims=(-maxabsz,maxabsz),seriescolor=:balance)
 xlabel!(plt,"v (m/s)")
 ylabel!(plt,"Observation Time (hr)")
 title!(plt,"Sun " * date_str * ":\nCCF-<CCF>+v⋅d<CCF>/dv")
 #xlims!(ccf_mid_velocity-RvSpectMLBase.max_bc,ccf_mid_velocity+RvSpectMLBase.max_bc)
 xlims!(plt,ccf_mid_velocity-2*line_width_50,ccf_mid_velocity+2*line_width_50)
 hline!(plt,[-Δt_fit,Δt_fit].*24, label=:none, color=:grey)
 hline!(plt,[0.], label=:none, color=:orange)
 ylims!(plt,-3.5,3.5)
 display(plt)
 savefig(joinpath(output_dir,"ccf_resid_heatmaps_sun_" * date_str * "_espresso_msf=" * string(round(msf,digits=1)) * ".png"))



#=
 all_spectra = nothing
 order_list_timeseries = nothing
 GC.gc()
end
=#

#EchelleCCFs.write_each_ccf_fits(map(s->s.metadata,all_spectra), v_grid, ccfs_espresso, ccf_vars_espresso,   fits_hdr=EchelleCCFs.make_ccf_fits_header(all_spectra[1].metadata),
#            line_list_filename=line_list_filename,total_vel=rvs_ccf_espresso, sigma_total_vel=σ_rvs_espresso )


#=
max_orders = 56:111
  order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=max_orders, recalc=true )

#linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",line_list_filename)
 line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity,
      Δv_to_avoid_tellurics = 1*RvSpectMLBase.max_bc+5*line_width_50+2*default_ccf_mask_v_width(EXPRES1D()),
      orders_to_use=max_orders, recalc=true, verbose=true)
=#
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

#=
outputs[""] =
outputs[""] =
outputs[""] =
outputs[""] =

=#

order_snr = mapreduce(obsid->map(ord->RvSpectMLBase.calc_snr(all_spectra[obsid],787:6214,ord)  ,1:size(first(all_spectra).flux,2)), hcat, 1:length(all_spectra))
outputs["order_snr"] = order_snr

save(outputs_filename, outputs)

#=

EchelleCCFs.write_each_ccf_fits(map(s->s.metadata,all_spectra), #= 100.0 .* EXPRESS USED cm/s =# v_grid, ccfs_espresso, ccf_vars_espresso,
                  fits_hdr=EchelleCCFs.make_ccf_fits_header(all_spectra[1].metadata),
                  line_list_filename=line_list_filename,# orders=161 .-good_orders,ccf_orders=order_ccfs,ccf_orders_var=order_ccf_vars,
                  total_vel=100.0 .* rvs_ccf_espresso, sigma_total_vel= 100.0 .* σ_rvs_espresso,
                  order_vels = 100.0 .* order_rvs_g, sigma_order_vels = 100.0 .* order_rv_std_g )
=#
|
