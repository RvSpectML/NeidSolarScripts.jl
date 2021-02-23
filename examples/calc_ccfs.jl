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


target_subdir = "good_days"   # USER: Replace with directory of your choice
  fits_target_str = "Sun"
  output_dir = "examples/output/"
  paths_to_search_for_param = [pwd(),"examples",joinpath(pkgdir(RvSpectMLBase),"..","RvSpectML","examples"), "/gpfs/scratch/jpn23/"]
  # NOTE: make_manifest does not update its paths_to_search when default_paths_to_search is defined here, so if you change the line above, you must also include "paths_to_search=default_paths_to_search" in the make_manifest() function call below
  pipeline_plan = PipelinePlan()
  dont_make_plot!(pipeline_plan, :movie)

reset_all_needs!(pipeline_plan)
#if need_to(pipeline_plan,:read_spectra)
if verbose println("# Finding what data files are avaliable.")  end
    eval(read_data_paths(paths_to_search=paths_to_search_for_param))
    @assert isdefined(Main,:neid_data_path)
    df_files = make_manifest(neid_data_path, target_subdir, NEID )

if verbose println("# Reading in customized parameters from param.jl.")  end
  idx_day_to_use = 1
    eval(code_to_include_param_jl(paths_to_search=paths_to_search_for_param))
    if match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1] ==  match(r"neidL1_(\d+)[T_](\d+)\.fits$", last(df_files_use.Filename))[1]
      date_str = match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1]
    else
      date_str = string(match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1]) * "-" * string(match(r"neidL1_(\d+)[T_](\d+)\.fits$", last(df_files_use.Filename))[1])
    end

#for idx_day_to_use in 1:size(df_files_solar_by_day,1)
   #df_files_use = df_files_solar_by_day[idx_day_to_use,:data] |> @orderby(_.bjd) |> @take(max_spectra_to_use) |> DataFrame

if verbose println("# Reading in ", size(df_files_use,1), " FITS files.")  end
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
    order_list_timeseries = extract_orders(all_spectra,pipeline_plan,  orders_to_use=good_orders,  recalc=true )

    line_list_filename = "G2.espresso.mas"
    linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",line_list_filename)

    line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity,
       Δv_to_avoid_tellurics = 1*max_bc+5*line_width_50_default+max_mask_scale_factor*default_ccf_mask_v_width(NEID1D()),#= orders_to_use=good_orders, =# recalc=true, verbose=true)

    CSV.write(joinpath("data","solar_line_list_espresso.csv"), line_list_espresso)
 end
 line_width_50 = line_width_50_default


#good_orders = orders_to_use_default(NEID2D())
good_orders = 55 .+ vcat(1:14, 17:24, 27:29, 37:38,40:41)
 order_list_timeseries = extract_orders(all_spectra,pipeline_plan,  orders_to_use=good_orders,  recalc=true )

 #msf = 1.4; fwtf = 1.5
 msf = 3.5; fwtf = 1.5
 println("# mask_scale_factor = ", msf, ", frac_width_to_fit =  ", fwtf, ". ")
 ((ccfs_espresso, ccf_vars_espresso), v_grid) = ccf_total(order_list_timeseries, line_list_espresso, pipeline_plan, #= mask_type=:gaussian,=#
   mask_scale_factor=msf, range_no_mask_change=5*line_width_50, ccf_mid_velocity=ccf_mid_velocity, v_step=200,
   v_max=5*line_width_50, calc_ccf_var=true, recalc=true)
 #fwtf = 2.0
 alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
 rvs_ccf_espresso = RvSpectML.calc_rvs_from_ccf_total(ccfs_espresso, ccf_vars_espresso, pipeline_plan, v_grid=v_grid, times = order_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
 σ_rvs_espresso = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))
 #scatter((order_list_timeseries.times.-first(order_list_timeseries.times)).*24, rvs_ccf_espresso, label=:none)
 times = order_list_timeseries.times.-first(order_list_timeseries.times)
 x = [times  ones(length(times)) ]
 linfit = (x'*x)\x'*(rvs_ccf_espresso)
 v_fit = x*linfit
 (times_binned, rvs_binned) = bin_times_and_rvs_max_Δt(times=times, rvs=rvs_ccf_espresso, Δt_threshold=5/(60*24))
 x = [times_binned  ones(length(times_binned)) ]
 v_fit_binned = x*linfit
 println("# RMS to linear fit: ", std(rvs_ccf_espresso.-v_fit), " m/s.  Binned: ", std(rvs_binned.-v_fit_binned))
 using Plots
 scatter( times.*24, rvs_ccf_espresso.-v_fit,  ms=1.5, label=:none)
 scatter!(times_binned.*24, rvs_binned.-v_fit_binned, ms=4.0, label=:none)
 xlabel!("Time (hr)")
 ylabel!("RV - linear fit (m/s)")
 title!("Sun " * date_str * " CCF ESPRESSO\nSlope = " * string(round(linfit[1]/24,digits=3))  * "m/s/hr" * " rms_5min = " * string(round(std(rvs_binned.-v_fit_binned),digits=3)) * "m/s")
 savefig("rvs_" * date_str * ".png")


ccfs_norm = EchelleCCFs.calc_normalized_ccfs( v_grid, ccfs_espresso, v_mid=ccf_mid_velocity, dv_min=2.5*line_width_50,dv_max= 2*RvSpectMLBase.max_bc)
 ccf_template = EchelleCCFs.calc_ccf_template(ccfs_norm, assume_normalized=true)
 plot(v_grid,ccfs_espresso, labels=:none)
 plot(v_grid,ccfs_norm , labels=:none)
 plot(v_grid,(ccfs_norm.-ccf_template)[:,range(1,size(ccfs_norm,2),step=2)], labels=:none)
 xlims!(ccf_mid_velocity-3*line_width_50,ccf_mid_velocity+3*line_width_50)
 heatmap(v_grid,1:size(ccfs_norm,2),(ccfs_norm.-ccf_template)')
 cols_to_fit = findlast(x->x<ccf_mid_velocity-3*line_width_50,v_grid):findfirst(x->x>ccf_mid_velocity+3*line_width_50,v_grid)
 heatmap(view(v_grid,cols_to_fit),1:size(ccfs_norm,2),view(ccfs_norm,cols_to_fit,:)'.-view(ccf_template,cols_to_fit,:)')
 xlabel!("v (m/s)")
 ylabel!("Observation ID")
 title!("Sun " * date_str * ":\nCCF - <CCF>, ESPRESSO mask + SNR weights")
 #xlims!(ccf_mid_velocity-RvSpectMLBase.max_bc,ccf_mid_velocity+RvSpectMLBase.max_bc)
 xlims!(ccf_mid_velocity-3*line_width_50,ccf_mid_velocity+3*line_width_50)
 savefig("ccf_heatmaps_sun_" * date_str * "_espresso_msf=" * string(round(msf,digits=1)) * ".png")

#=
 all_spectra = nothing
 order_list_timeseries = nothing
 GC.gc()
end
=#

#EchelleCCFs.write_each_ccf_fits(map(s->s.metadata,all_spectra), v_grid, ccfs_espresso, ccf_vars_espresso,   fits_hdr=EchelleCCFs.make_ccf_fits_header(all_spectra[1].metadata),
#            line_list_filename=line_list_filename,total_vel=rvs_ccf_espresso, sigma_total_vel=σ_rvs_espresso )



#=
max_orders = 12:83
  order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=max_orders, recalc=true )

linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",line_list_filename)
 line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity,
      Δv_to_avoid_tellurics = 1*RvSpectMLBase.max_bc+5*line_width_50+2*default_ccf_mask_v_width(EXPRES1D()),
      orders_to_use=max_orders, recalc=true, verbose=true)
=#

good_orders = orders_to_use_default(NEID2D())

(order_ccfs, order_ccf_vars, v_grid_order_ccfs) = ccf_orders(order_list_timeseries, line_list_espresso, pipeline_plan, mask_scale_factor=msf, ccf_mid_velocity=ccf_mid_velocity,
 v_step=100, v_max=5*line_width_50, orders_to_use=good_orders, calc_ccf_var=true,recalc=true)

obs_id = 1
heatmap(v_grid_order_ccfs, good_orders, (order_ccfs[:,:,obs_id]./mean(order_ccfs[:,:,obs_id],dims=1))')
xlabel!("v (m/s)")
ylabel!("Order index")
title!("Sun " * date_str * ":\nCCF ESPRESSO mask + SNR weights")
#xlims!(ccf_mid_velocity-RvSpectMLBase.max_bc,ccf_mid_velocity+RvSpectMLBase.max_bc)
xlims!(ccf_mid_velocity-3*line_width_50,ccf_mid_velocity+3*line_width_50)
savefig("order_ccf_heatmaps_sun_" * date_str * "_espresso_msf=" * string(round(msf,digits=1)) * ".png")


alg_fit_rv2 =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
order_rvs_g = zeros(length(all_spectra), length(good_orders))
order_rv_std_g = zeros(length(all_spectra), length(good_orders))
order_rvs_t = zeros(length(all_spectra), length(good_orders))
order_rv_std_t = zeros(length(all_spectra), length(good_orders))
for (i,ord) in enumerate(good_orders)
   println("# Order: ", i)
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

#=

EchelleCCFs.write_each_ccf_fits(map(s->s.metadata,all_spectra), #= 100.0 .* EXPRESS USED cm/s =# v_grid, ccfs_espresso, ccf_vars_espresso,
                  fits_hdr=EchelleCCFs.make_ccf_fits_header(all_spectra[1].metadata),
                  line_list_filename=line_list_filename,# orders=161 .-good_orders,ccf_orders=order_ccfs,ccf_orders_var=order_ccf_vars,
                  total_vel=100.0 .* rvs_ccf_espresso, sigma_total_vel= 100.0 .* σ_rvs_espresso,
                  order_vels = 100.0 .* order_rvs_g, sigma_order_vels = 100.0 .* order_rv_std_g )
=#
|
