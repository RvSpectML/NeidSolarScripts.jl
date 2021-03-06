verbose = true
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectMLBase
 using EchelleInstruments, EchelleInstruments.NEID
 using RvSpectML
 if verbose println("# Loading NeidSolarScripts")    end
 using NeidSolarScripts
 #using NeidSolarScripts.SolarRotation
 #using NeidSolarScripts.DifferentialExtinction
 if verbose   println("# Loading other packages")    end
 using CSV, DataFrames, Query, Dates
 using JLD2, FileIO
 #using StatsBase, Statistics, Dates

fits_target_str = "Sun"
  #target_subdir = "/mnt/data_simons/NEID/DRPv0.7-fixedflatfielding2"   # USER: Replace with directory of your choice
  #output_dir = "output/"
  target_subdir = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/outputs/20210330"
  output_dir = ""
 if !@isdefined(target_subdir)
      target_subdir = ""
 end
 if !@isdefined(create_missing_continuum_files)
     create_missing_continuum_files = false
 end
 paths_to_search_for_param = [pwd(),joinpath(dirname(pathof(NeidSolarScripts)),"..")]

 if verbose println("# Finding what data files are avaliable.")  end
 eval(read_data_paths(paths_to_search=paths_to_search_for_param))
 @assert isdefined(Main,:neid_data_path)
 @assert isdefined(Main,:output_dir)
 output_dir = ""
 output_path = joinpath(neid_data_path, target_subdir,output_dir)
 manifest_filename = joinpath(output_path,"manifest.csv")
 manifest_calib_filename = joinpath(output_path,"manifest_calib.csv")

if isfile(manifest_filename)
    df_files  = CSV.read(manifest_filename, DataFrame)
    @assert size(df_files,1) >= 1
    @assert hasproperty(df_files,:Filename)
    @assert hasproperty(df_files,:target)
    @assert hasproperty(df_files,:bjd)
    @assert hasproperty(df_files,:ssbz)
    @assert hasproperty(df_files,:exptime)
    @assert hasproperty(df_files,:alt_sun)
    @assert hasproperty(df_files,:Δfwhm²)
    @assert hasproperty(df_files,:order_snrs)
    if eltype(df_files[!,:order_snrs]) == String
        df_files[!,:order_snrs] = map(i->parse.(Float64,split(df_files[i,:order_snrs][2:end-1],',')),1:size(df_files,1))
    end
    @assert eltype(df_files[!,:order_snrs]) == Vector{Float64}
else
    @error "Need to make manifest at $manifest_filename"
    #df_files = NEID.make_manifest(joinpath(neid_data_path, target_subdir))
    #CSV.write(manifest_filename, df_files)
end


paths_to_search_for_param = [pwd(),"..","examples",joinpath(pkgdir(RvSpectMLBase),"..","RvSpectML","examples"), "/gpfs/scratch/jpn23/"]
  # NOTE: make_manifest does not update its paths_to_search when default_paths_to_search is defined here, so if you change the line above, you must also include "paths_to_search=default_paths_to_search" in the make_manifest() function call below
  #outputs = Dict{String,Any}()
  pipeline_plan = PipelinePlan()
  dont_make_plot!(pipeline_plan, :movie)
  if verbose println("# Reading in customized parameters from param.jl.")  end
   if !@isdefined(idx_day_to_use)
       idx_day_to_use = 1
   end
   eval(code_to_include_param_jl(paths_to_search=paths_to_search_for_param))
   #outputs["df_files_use"] = df_files_use

if match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1] ==  match(r"neidL1_(\d+)[T_](\d+)\.fits$", last(df_files_use.Filename))[1]
      date_str = match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1]
    else
      date_str = string(match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1]) * "-" * string(match(r"neidL1_(\d+)[T_](\d+)\.fits$", last(df_files_use.Filename))[1])
   end

#=
   outputs_filename = joinpath(output_dir,"solar_" * date_str * "_new.jld2")
   if isfile(outputs_filename) && false
     times_already_processed = load(outputs_filename, "times")
     files_in_day_to_process = size(df_files_solar_by_day.data[idx_day_to_use],1)
      if files_in_day_to_process == length(times_already_processed)
         println("# Already processed all ", length(times_already_processed), " files for ", date_str)
         exit()
      end
   end
=#

#for idx_day_to_use in 1:size(df_files_solar_by_day,1)
   #df_files_use = df_files_solar_by_day[idx_day_to_use,:data] |> @orderby(_.bjd) |> @take(max_spectra_to_use) |> DataFrame

#=
if verbose println("# Reading in ", size(df_files_use,1), " FITS files for ", date_str, ".")  end
    @time all_spectra = map(row->NEID.read_data(row), eachrow(df_files_use))
    GC.gc()
    dont_need_to!(pipeline_plan,:read_spectra)
=#

all_spectra = Spectra2DBasic{Float64, Float32, Float32, Matrix{Float64}, Matrix{Float32}, Matrix{Float32}, NEID2D}[]
for (i,row) in enumerate(eachrow(df_files_use))
    m = match(r"(neidL1_\d+[T_]\d+)\.fits$", row.Filename)
    #continuum_filename = joinpath(neid_data_path,target_subdir,"output","continuum", m[1] * "_continuum=afs.jld2")
    continuum_filename = joinpath(neid_data_path,target_subdir,output_dir,"continuum", m[1] * "_continuum=afs.jld2")
    if !(isfile(continuum_filename)||islink(continuum_filename))
        println("# Couldn't find matching continuum file", continuum_filename)
        continue
    end
    spec = NEID.read_data(row)
    continuum = load(continuum_filename, "continuum")
    spec.flux ./= continuum
    spec.var ./= continuum.^2
    push!(all_spectra,spec)
    #if i>=20  break   end
end
GC.gc()
dont_need_to!(pipeline_plan,:read_spectra)


line_width_50_default = 7.9e3
 lsf_width = 3.0e3
 max_mask_scale_factor = max(lsf_width,line_width_50_default/sqrt(8*log(2)))/default_ccf_mask_v_width(NEID2D()) # 4.0
 max_bc = RvSpectMLBase.max_bc
 #max_bc = RvSpectMLBase.max_bc_earth_rotation
 line_list_filename = joinpath("data","solar_line_list_espresso.csv")
 max_orders = min_order(NEID2D()):max_order(NEID2D())
 good_orders = orders_to_use_default(NEID2D())
 orders_to_use = max_orders
 if isfile(line_list_filename)
   println("# Reading ", line_list_filename)
   line_list_espresso = CSV.read(line_list_filename, DataFrame)
   dont_need_to!(pipeline_plan,:clean_line_list_tellurics)
 else
    println("# Can't find ", line_list_filename, ".  Trying ESPRESSO line list.")
    orders_to_use = good_orders
    order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=orders_to_use, recalc=true )

    line_list_filename = "G2.espresso.mas"
    linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",line_list_filename)

    line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity,
       Δv_to_avoid_tellurics = 2*max_bc+5*line_width_50_default+max_mask_scale_factor*default_ccf_mask_v_width(NEID2D()), orders_to_use=orders_to_use, recalc=true, verbose=true)

    #CSV.write(custom_line_list_filename, line_list_espresso)
 end
 #outputs["line_list_espresso"] = line_list_espresso

line_width_50 = line_width_50_default
 maxΔfwhm² = -0.569375
 @assert maximum(df_files_use.Δfwhm²) < maxΔfwhm²
 Δfwhm = 1000.0 .* sqrt.(maxΔfwhm².-df_files_use.Δfwhm²[1:length(all_spectra)])  # How much to increase fwhm by to acheive uniform fwhm


msf = lsf_width/default_ccf_mask_v_width(NEID2D()); fwtf = 0.5  # using LSF width
 order_list_timeseries = extract_orders(all_spectra,pipeline_plan,  orders_to_use=orders_to_use,  remove_bad_chunks=false, recalc=true )

@time (order_ccfs, order_ccf_vars, v_grid_order_ccfs) = ccf_orders(order_list_timeseries, line_list_espresso, pipeline_plan,
    mask_type=:gaussian, Δfwhm=Δfwhm,
    mask_scale_factor=msf, range_no_mask_change=5*line_width_50, ccf_mid_velocity=ccf_mid_velocity, v_step=155,
    v_max=max(5*line_width_50,2*max_bc), orders_to_use=orders_to_use, allow_nans=true, calc_ccf_var=true,
    recalc=true)

orders_to_use2 = orders_to_use[map(i->!iszero(order_ccfs[:,i,:]),1:length(orders_to_use))]
order_list_timeseries2 = extract_orders(all_spectra,pipeline_plan,  orders_to_use=orders_to_use2,  remove_bad_chunks=false, recalc=true )

ccf_dir = joinpath(neid_data_path,target_subdir,"ccfs")
if !isdir(ccf_dir)
   mkdir(ccf_dir)
end

#daily_ccf_filename = joinpath(neid_data_path,target_subdir,"output","ccfs", date_str * "_ccfs=default.jld2")
daily_ccf_filename = joinpath(neid_data_path,target_subdir,"ccfs", date_str * "_ccfs=default.jld2")
  jldopen(daily_ccf_filename, "w") do f
    f["v_grid"] = v_grid_order_ccfs
    f["order_ccfs"] = order_ccfs
    f["order_ccf_vars"] = order_ccf_vars
    f["orders_to_use"] = orders_to_use
    f["manifest"] = df_files_use
    f["daily_ccf_filename"] = daily_ccf_filename
  end


for (i,row) in enumerate(eachrow(df_files_use))
    m = match(r"(neidL1_\d+[T_]\d+)\.fits$", row.Filename)
    #ccf_filename = joinpath(neid_data_path,target_subdir,"output","ccfs", m[1] * "_ccfs=default.jld2")
    ccf_filename = joinpath(neid_data_path,target_subdir,"ccfs", m[1] * "_ccfs=default.jld2")
    jldopen(ccf_filename, "w") do f
        f["v_grid"] = v_grid_order_ccfs
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
  outputs["v_grid"] = v_grid
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
