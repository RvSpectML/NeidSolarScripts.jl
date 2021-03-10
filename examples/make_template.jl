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

alt_sun = map(i->calc_solar_alt(df_files_use[i,:bjd]),1:size(df_files_use,1))
obs_idx_max_alt_sun = argmax(alt_sun)
time_max_alt_sun = df_files_use[obs_idx_max_alt_sun,:bjd]
df_files_use.Δt = df_files_use[!,:bjd].-time_max_alt_sun
max_Δt = 0.5/24
df_files_use = df_files_use |> @filter( abs(_.Δt) <= max_Δt ) |> DataFrame

#for idx_day_to_use in 1:size(df_files_solar_by_day,1)
   #df_files_use = df_files_solar_by_day[idx_day_to_use,:data] |> @orderby(_.bjd) |> @take(max_spectra_to_use) |> DataFrame

if verbose println("# Reading in ", size(df_files_use,1), " FITS files for ", date_str, ".")  end
    @time all_spectra = map(row->NEID.read_data(row), eachrow(df_files_use))
    GC.gc()
    dont_need_to!(pipeline_plan,:read_spectra)

good_orders = 56:111 # Everything with λ and flux
 # good_orders = orders_to_use_default(NEID2D())
 #order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=max_orders, recalc=true )
 #order_list_timeseries = extract_orders(all_spectra,pipeline_plan,  orders_to_use=good_orders,  remove_bad_chunks=false, recalc=true )

df_orders_pixels = NEID.make_good_orders_pixels_df(all_spectra ; orders=good_orders)
chunk_list_timeseries = RvSpectMLBase.make_chunk_list_timeseries_from_orders_pixels_df(all_spectra,NEID2D(), df_orders_pixels)



need_to!(pipeline_plan, :template)
 if need_to(pipeline_plan, :template)  # Compute order CCF's & measure RVs
    if verbose println("# Making template spectra.")  end
    #@assert !need_to(pipeline_plan,:extract_orders)
    map(i->chunk_list_timeseries.metadata[i][:rv_est] = 0.0, 1:length(chunk_list_timeseries) )
    GC.gc()   # run garbage collector for deallocated memory
    #map(i->order_list_timeseries.metadata[i][:rv_est] = rvs_ccf[i]-mean(rvs_ccf), 1:length(order_list_timeseries) )
    # Smothing is broken with GP interpolation.  Need to fix.  In mean time, here's a Sinc interpolation workaround
    @time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(chunk_list_timeseries, smooth_factor=2.0)
    #@time ( spectral_orders_matrix, f_mean, var_mean, deriv, deriv2, order_grids )  = RvSpectML.make_template_spectra(order_list_timeseries, alg=:Sinc)
    if save_data(pipeline_plan, :template)
       using JLD2, FileIO
       save(joinpath(output_dir,target_subdir * "_matrix.jld2"), Dict("λ"=>spectral_orders_matrix.λ,"spectra"=>spectral_orders_matrix.flux,"var_spectra"=>spectral_orders_matrix.var,"flux_template"=>f_mean,"var"=>var_mean, "dfluxdlnλ_template"=>deriv,"d²fluxdlnλ²_template"=>deriv2))
    end
    dont_need_to!(pipeline_plan, :template);
 end

#=
chunkid = 21
 idx = spectral_orders_matrix.chunk_map[chunkid]
 plt = plot(spectral_orders_matrix.λ[idx],f_mean[idx],markersize=1.0,label="Template")
 plot!(plt,chunk_list_timeseries[1][chunkid].λ,chunk_list_timeseries[1][chunkid].flux, markersize=1.2,label="Obs  1")
 plot!(plt,chunk_list_timeseries[10][chunkid].λ,chunk_list_timeseries[10][chunkid].flux, markersize=1.2,label="Obs 10")
 plot!(plt,chunk_list_timeseries[20][chunkid].λ,chunk_list_timeseries[20][chunkid].flux, markersize=1.2,label="Obs 20")
 #plot!(plt,spectral_orders_matrix.λ[idx],deriv[idx]./mean(abs.(deriv[idx]))/2,markersize=1.1,label="Deriv")
 #plot!(plt,spectral_orders_matrix.λ[idx],deriv2[idx]./mean(abs.(deriv2[idx]))/2,markersize=1.1,label="Deriv 2")
 xlabel!("λ (Å)")
 ylabel!("f(λ)")
 title!("Template spectrum for chunk " * string(chunkid) )
 xlims!(spectral_orders_matrix.λ[idx[floor(Int,1+0.2*length(idx))]],spectral_orders_matrix.λ[idx[floor(Int,1+0.35*length(idx))]])

chunkid = 21
 idx = spectral_orders_matrix.chunk_map[chunkid]
 plt = plot(spectral_orders_matrix.λ[idx],(f_mean[idx].-1.0)./mean(abs.((f_mean[idx].-1.0))),markersize=1.0,label="Template")
 plot!(plt,spectral_orders_matrix.λ[idx],deriv[idx]./mean(abs.(deriv[idx])),markersize=1.1,label="Deriv")
 plot!(plt,spectral_orders_matrix.λ[idx],deriv2[idx]./mean(abs.(deriv2[idx])),markersize=1.1,label="Deriv 2")
 xlabel!("λ (Å)")
 ylabel!("f(λ), f'(λ), f''(λ), all standardized")
 title!("Template spectrum for chunk " * string(chunkid) )
 xlims!(spectral_orders_matrix.λ[idx[floor(Int,1+0.2*length(idx))]],spectral_orders_matrix.λ[idx[floor(Int,1+0.35*length(idx))]])
=#

line_width_50_default = 7.9e3
  lsf_width = 3.0e3
  max_mask_scale_factor = 4.0
  max_bc = RvSpectMLBase.max_bc
  #max_bc = RvSpectMLBase.max_bc_earth_rotation
  #if isfile(joinpath("data","solar_line_list_espresso.csv"))
   # line_list_espresso = CSV.read(joinpath("data","solar_line_list_espresso.csv"), DataFrame)
   # dont_need_to!(pipeline_plan,:clean_line_list_tellurics)
  #else
     #max_orders = min_order(NEID2D()):max_order(NEID2D())
     #good_orders = orders_to_use_default(NEID2D())
     #order_list_timeseries = extract_orders(all_spectra,pipeline_plan, orders_to_use=max_orders, recalc=true )
     #order_list_timeseries = extract_orders(all_spectra,pipeline_plan,  orders_to_use=good_orders,  remove_bad_chunks=false, recalc=true )

     line_list_filename = "G2.espresso.mas"
     linelist_for_ccf_fn_w_path = joinpath(pkgdir(EchelleCCFs),"data","masks",line_list_filename)

     line_list_espresso = prepare_line_list(linelist_for_ccf_fn_w_path, all_spectra, pipeline_plan, v_center_to_avoid_tellurics=ccf_mid_velocity,
        Δv_to_avoid_tellurics = 2*max_bc+5*line_width_50_default+max_mask_scale_factor*default_ccf_mask_v_width(NEID1D()), orders_to_use=good_orders, recalc=true, verbose=true)

     #CSV.write(joinpath("data","solar_line_list_espresso.csv"), line_list_espresso)
  #end
  #outputs["line_list_espresso"] = line_list_espresso
  line_width_50 = line_width_50_default


alt_sun = map(i->calc_solar_alt(all_spectra[i].metadata[:bjd]),1:length(all_spectra))
  obs_idx_max_alt_sun = argmax(alt_sun)
  time_max_alt_sun = chunk_list_timeseries.times[obs_idx_max_alt_sun]
  order_snrs = [RvSpectMLBase.calc_snr(chunk_list_timeseries[obsid].data[c]) for c in 1:length(chunk_list_timeseries[1].data),
                     obsid in 1:length(chunk_list_timeseries) ]
  order_weights =  (order_snrs[:,obs_idx_max_alt_sun]./order_snrs).^2
  for obsid in 1:length(chunk_list_timeseries)
     for ch in 1:length(chunk_list_timeseries[1].data)
         chunk_list_timeseries[obsid].data[ch].flux .*= order_weights[ch,obsid] #./mean_order_mean_flux[ch]
         chunk_list_timeseries[obsid].data[ch].var .*= (order_weights[ch,obsid]).^2 # ./mean_order_mean_flux[ch]).^2
     end
  end

msf = 0.62; fwtf = 1.5  # one pixel because tophat not using NEID's default v width
 #msf = 1.24; fwtf = 1.5   # two pixels
  println("# mask_scale_factor = ", msf, ", frac_width_to_fit =  ", fwtf, ". ")
  dont_need_to!(pipeline_plan,:extract_orders)
  ((ccfs_espresso, ccf_vars_espresso), v_grid) = ccf_total(chunk_list_timeseries, line_list_espresso, pipeline_plan,
    #mask_type=:gaussian,
    mask_scale_factor=msf, range_no_mask_change=5*line_width_50, ccf_mid_velocity=ccf_mid_velocity, v_step=155,
    v_max=max(5*line_width_50,2*max_bc), allow_nans=true, calc_ccf_var=true, recalc=true)

fwtf = 1.5; println("# mask_scale_factor = ", msf, ", frac_width_to_fit =  ", fwtf, ". ")
  alg_fit_rv =  EchelleCCFs.MeasureRvFromCCFGaussian(frac_of_width_to_fit=fwtf)
  rvs_ccf_espresso = RvSpectML.calc_rvs_from_ccf_total(ccfs_espresso, ccf_vars_espresso, pipeline_plan, v_grid=v_grid, times = chunk_list_timeseries.times, recalc=true, bin_nightly=true, alg_fit_rv=alg_fit_rv)
  σ_rvs_espresso = deepcopy(read_cache(pipeline_plan,:σ_rvs_ccf_total))



need_to!(pipeline_plan,:fit_lines)
 if need_to(pipeline_plan,:fit_lines)
    if verbose println("# Performing fresh search for lines in template spectra.")  end
    cl = ChunkList(map(grid->ChunkOfSpectrum(spectral_orders_matrix.λ,f_mean, var_mean, grid), spectral_orders_matrix.chunk_map), ones(Int64,length(spectral_orders_matrix.chunk_map)))
    # We're done with the spectral_orders_matrix, so we can release the memory now
    #spectral_orders_matrix = nothing
    #GC.gc()
    need_to!(pipeline_plan,:template)
    #lines_in_template_logy = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(line_width=line_width_50,min_deriv2=0.5, use_logλ=true, use_logflux=true), verbose=true)  # TODO: Automate threshold for finding a line
    #line_width_50_default = 7.9e3
    #line_width_50 = line_width_50_default
    line_finder_plan = RvSpectML.LineFinder.LineFinderPlan(line_width=line_width_50,min_deriv2=0.5, use_logλ=true, use_logflux=false)
    lines_in_template = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=line_finder_plan, verbose=true)  # TODO: Automate threshold for finding a line

    if verbose println("# Finding above lines in all spectra.")  end
    rvs_ccf = rvs_ccf_espresso
    @time fits_to_lines = RvSpectML.LineFinder.fit_all_lines_in_chunklist_timeseries(chunk_list_timeseries, lines_in_template, plan=line_finder_plan, rvs=rvs_ccf )

    #if save_data(pipeline_plan,:fit_lines)
       using CSV
       CSV.write(joinpath(output_dir, date_str * "_linefinder_lines.csv"), lines_in_template )
       CSV.write(joinpath(output_dir, date_str * "_linefinder_line_fits.csv"), fits_to_lines )
       #CSV.write(joinpath(output_dir,target_subdir * "_linefinder_line_fits_clean.csv"), lines_to_try )
    #end
    dont_need_to!(pipeline_plan,:fit_lines);
 end

function calc_lin_fit_coeff(times::AA1, y::AA2)  where { T1<:Real, T2<:Real, AA1<:AbstractArray{T1,1}, AA2<:AbstractArray{T2,1} }
    x = [ ones(length(times)) times ]
    linfit = (x'*x)\x'*y
 end

function calc_quad_fit_coeff(times::AA1, y::AA2)  where { T1<:Real, T2<:Real, AA1<:AbstractArray{T1,1}, AA2<:AbstractArray{T2,1} }
   x = [ ones(length(times)) times  times.^2 ]
   quadfit = (x'*x)\x'*y
end

function calc_ddt_coeff(times::AA1, y::AA2)  where { T1<:Real, T2<:Real, AA1<:AbstractArray{T1,1}, AA2<:AbstractArray{T2,1} }
   linfit = calc_lin_fit_coeff(times,y)
   slope = linfit[2]
end

using Printf

function select_line_fits_with_good_depth_width_slope(line_fits_df::DataFrame, quantile_threshold::Real; verbose::Bool = false, write_csv::Bool = false )
    conv_frac = line_fits_df |> @groupby(_.line_id) |> @map( { line_id=first(_.line_id), frac_converged=mean(_.fit_converged) } ) |> DataFrame

    fit_distrib = line_fits_df |> @filter(_.fit_converged == 1) |> @groupby(_.line_id) |>
             @map( { median_a=median(_.fit_a), median_b=median(_.fit_b), median_depth=median(_.fit_depth), median_σ²=median(_.fit_σ²), median_λc=median(_.fit_λc),
                    std_a=std(_.fit_a), std_b=std(_.fit_b), std_depth=std(_.fit_depth), std_σ²=std(_.fit_σ²), std_λc=std(_.fit_λc),
#                    Δt = df_files_use.Δt[_.obs_idx],
#                    ddt_depth=calc_ddt_coeff(df_files_use.Δt,_.fit_depth),
#                    ddt_σ²=calc_ddt_coeff(Δt,_.fit_σ²),
#                    ddt_λc=calc_ddt_coeff(Δt,_.fit_λc),
                    line_id=first(_.line_id) }) |>
             @join(conv_frac, _.line_id, _.line_id, {line_id=_.line_id, median_a=_.median_a, median_b=_.median_b, median_depth=_.median_depth, median_σ²=_.median_σ², median_λc=_.median_λc,
                    std_a=_.std_a, std_b=_.std_b, std_depth=_.std_depth, std_σ²=_.std_σ², std_λc=_.std_λc,
#                    ddt_depth=_.ddt_depth, ddt_σ²=_.ddt_σ², ddt_λc = _.ddt_λc,
                    frac_converged=__.frac_converged } ) |>
             @filter(_.frac_converged > 0.5 ) |> DataFrame

    std_depth_treshold = quantile(fit_distrib.std_depth,quantile_threshold)
    #median_σ_width_treshold = quantile(sqrt.(fit_distrib.median_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,1-quantile_threshold)
    median_σ_width_treshold_lo = 2000
    #median_σ_width_treshold_lo = 100
    #median_σ_width_treshold = line_width_50/3
    std_σ_width_treshold = quantile(sqrt.(fit_distrib.std_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,quantile_threshold)
    #std_depth_treshold = quantile(fit_distrib.std_depth,quantile_threshold)
    #median_σ_width_treshold_lo = quantile(sqrt.(fit_distrib.median_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,(1-quantile_threshold)/2)
    median_σ_width_treshold_hi = quantile(sqrt.(fit_distrib.median_σ²)./fit_distrib.median_λc.*RvSpectML.speed_of_light_mps,1-(1-quantile_threshold)/2)
    std_b_treshold = quantile(fit_distrib.std_b,quantile_threshold)
    std_a_treshold = quantile(fit_distrib.std_a,quantile_threshold)
    good_lines_alt = fit_distrib |>
       @filter( 0.05 <= _.median_depth <= 1.0 ) |>
       #@filter( sqrt.(_.median_σ²)./_.median_λc.*RvSpectML.speed_of_light_mps >= median_σ_width_treshold_lo ) |>
       @filter( sqrt.(_.std_σ²)./_.median_λc.*RvSpectML.speed_of_light_mps <= std_σ_width_treshold ) |>
       #@filter(  median_σ_width_treshold_lo <= sqrt.(_.median_σ²)./_.median_λc.*RvSpectML.speed_of_light_mps <= median_σ_width_treshold_hi ) |>
       @filter( _.std_depth <= std_depth_treshold ) |>
       @filter( _.std_b < std_b_treshold) |>
       DataFrame
    if verbose
       println("# Found ", size(good_lines_alt,1), " good lines (std_depth_width_slope), rejected ", size(fit_distrib,1)-size(good_lines_alt,1), " lines.")
    end
    if write_csv
       val_str = Printf.@sprintf("%1.2f",quantile_threshold)
       CSV.write(joinpath(output_dir,date_str * "_good_lines_fit_quant=" * val_str * ".csv"), good_lines_alt )
    end
    return good_lines_alt
 end

function calc_rvs_from_line_by_line_fits_alt(fits_to_lines::DataFrame, threshold_alt::Real)
    good_lines_alt = select_line_fits_with_good_depth_width_slope(fits_to_lines,threshold_alt,write_csv=true)
    df2 = fits_to_lines |> @join(good_lines_alt, _.line_id, _.line_id, {obs_id=_.obs_idx, line_id=_.line_id, Δv=(_.fit_λc.-__.median_λc).*RvSpectML.speed_of_light_mps./__.median_λc, weight=(min(__.median_λc/__.std_λc, RvSpectML.speed_of_light_mps/750))^2
    #= (__.median_λc/(RvSpectML.speed_of_light_mps*__.std_λc))^2 =# }) |>
             @groupby(_.obs_id) |> @map({obs_id=first(_.obs_id), Δv= sum(_.Δv.*_.weight)./sum(_.weight) }) |> DataFrame
    rms_alt = std(df2.Δv)
 end


good_lines_alt = select_line_fits_with_good_depth_width_slope(fits_to_lines,0.8) #0.8)
df2 = fits_to_lines |> @join(good_lines_alt, _.line_id, _.line_id, {obs_id=_.obs_idx, line_id=_.line_id, Δv=(_.fit_λc.-__.median_λc).*RvSpectML.speed_of_light_mps./__.median_λc, weight=(min(__.median_λc/__.std_λc, RvSpectML.speed_of_light_mps/750))^2 }) |>
          @groupby(_.obs_id) |> @map({obs_id=first(_.obs_id), Δv= sum(_.Δv.*_.weight)./sum(_.weight) }) |> DataFrame

using Plots
#thresholds = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
thresholds = [ 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0 ]
 #thresholds = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 ]
 #thresholds = [0.75, 0.8, 0.9, 0.95, 1.0 ]
  plt = plot()
  #out_vec_calc_rvs_from_line_by_line_fits_stdv = map(t->calc_rvs_from_line_by_line_fits_stdv(fits_to_lines,t),thresholds)
  #scatter!(plt,thresholds,out_vec_calc_rvs_from_line_by_line_fits_stdv,color=1,label="σ_λ")
  #plot!(plt,thresholds,out_vec_calc_rvs_from_line_by_line_fits_stdv,color=1,label=:none)
  out_vec_calc_rvs_from_line_by_line_fits_alt = map(t->calc_rvs_from_line_by_line_fits_alt(fits_to_lines,t),thresholds)
  scatter!(plt,thresholds,out_vec_calc_rvs_from_line_by_line_fits_alt,color=2,label="σ_other")
  plot!(plt,thresholds,out_vec_calc_rvs_from_line_by_line_fits_alt,color=2,label=:none)
  xlabel!("Threshold for accepting lines")
  ylabel!("RMS v (m/s)")
  savefig("rms_rv_vs_line_accept_threshold" * date_str * ".png")
  #display(plt)

#hpf_mask_data = readdlm(joinpath(pwd(),"data","HPF_TelMask_interpolated.txt"))
