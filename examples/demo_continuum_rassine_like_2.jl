using NeidSolarScripts    # For continuum, solar rotaiton, differential rotation,...
 using RvSpectMLBase      # For spectra data structures
 using EchelleInstruments # For reading spectra
 using JLD2, FileIO       # For outputs
 using ArgParse           # For command line arguments

 # using CSV, DataFrames, Query, Dates, SHA
 #=
 println("# Parsing arguments...")
    function parse_commandline()
      s = ArgParseSettings( description = "Calculate continuum normalization from NEID L1 FITS files.")
      add_arg_group!(s, "Files", :argg_files)
      @add_arg_table! s begin
          "manifest"
              help = "Manifest file (CVS) containing FITS files to analyze.\nExpects columns named Filename, bjd, and target."
              arg_type = String
              default = "manifest.csv"
              #required = true
          "output"
              help = "Filename for output CCFs (jld2)"
              arg_type = String
              default = "daily_ccfs.jld2"
          "--overwrite"
             help = "Specify it's ok to overwrite the output file."
             action = :store_true
       end
       add_arg_group!(s, "Normalization parameters", :argg_norm_param)
       @add_arg_table! s begin
         "--orders_to_use"
            help = "First and last order _index_ to compute CCFs for"
            nargs = 2
            arg_type = Int64
            #default = [ min_order(NEID2D()), max_order(NEID2D()) ]
            #default = [ first(orders_to_use_default(NEID2D())), last(orders_to_use_default(NEID2D())) ]
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
            default = 6.0
          "--max_airmass"
             help = "Maximum airmass."
             arg_type = Float64
             #default = 10.0
          "--min_snr_factor"
             help = "Minimum SNR relative to max_snr."
             arg_type = Float64
             #default = 0.0
          "--max_snr"
             help = "Specify max_snr manually."
             arg_type = Float64
          "--start_time"
             help = "Specify daily start time for CCFs to be processed (HH MM)"
             nargs = 2
             arg_type = Int64
             #default = [ min_order(NEID2D()), max_order(NEID2D()) ]
             default = [0, 0]
         "--stop_time"
            help = "Specify daily stop time for CCFs to be processed (HH MM)"
            nargs = 2
            arg_type = Int64
            #default = [ min_order(NEID2D()), max_order(NEID2D()) ]
            default = [23, 59]
          "--max_num_spectra"
             help = "Maximum number of spectra to process."
             arg_type = Int
             default = 250
       end

      return parse_args(s)
  end
  args = parse_commandline()
  =#

sed_path = "/home/eford/Code/RvSpectMLEcoSystem/NeidSolarScripts/data"
   sed = Continuum.read_master_sed_neid(path=sed_path)

fits_path = "/mnt/data_simons/NEID/DRPv0.7-fixedflatfielding2/"
   fits_fn = "neidL1_20210104T182500.fits"
   continuum_afs_path = joinpath(fits_path,"output","continuum")

   continuum_afs_fn = "neidL1_20210104T182500_continuum=afs.jld2"
   continuum_ras_path = joinpath(fits_path,"output","continuum")
   continuum_ras_fn = "neidL1_20210104T182500_continuum=ras1.jld2"
   write_output = false
   spec = NEID.read_data(joinpath(fits_path,fits_fn))

min_order_for_continuum = min_order(NEID2D())
  max_order_for_continuum = max_order(NEID2D())
  orders_to_use = min_order_for_continuum:max_order_for_continuum
  (f_norm, var_norm) = Continuum.normalize_by_sed(spec.λ,spec.flux,spec.var, sed; orders_to_use=orders_to_use)

using Plots
order_idx = 110
 plt = scatter(spec.λ[:,order_idx],f_norm[:,order_idx],ms=1.5,legend=:none)
 #plot!(plt, lambda,poly, color=:blue, lw=2)
 #plt = scatter(lambda,f_obs./(sed_view.*coeff.(lambda)),ms=1.5,legend=:none)
 #plt
#savefig("fig.png")
 ylims!(0.,2.0)
true


 function get_pixel_range_for_continuum(λ::AA1, ord_idx::Integer) where { T1<:Real, AA1<:AbstractArray{T1,2} }
   inst = get_inst(spec)
   pix = min_pixel_in_order(get_inst(spec)):max_pixel_in_order(get_inst(spec))
 end
 #f_norm = spec.flux ./sed
 #xlabel!("λ")
 #ylabel!("Flux/SED")
 #savefig("fig.png")

#f_norm = spec.flux
@time (anchors, continuum_output, f_filtered) = Continuum.calc_continuum(spec.λ,f_norm,var_norm,stretch_factor=5.0, ν=1.0,
               merging_threshold=0.25, orders_to_use=orders_to_use, get_pixel_range=get_pixel_range_for_continuum,
               verbose=false)


order_idx = 90
  lambda_obs = view(spec.λ,:,order_idx)
  f_obs = view(f_norm,:,order_idx)
  var_obs = view(var_norm,:,order_idx)
  (anchors, continuum, f_filtered_internal) = Continuum.calc_continuum(lambda_obs,f_obs,var_obs,
       stretch_factor=5.0, ν=1.0, merging_threshold=0.25,  min_R_factor = 100.0, smoothing_half_width = 40, verbose=false)
  f_filtered = f_filtered_internal
  scatter(lambda_obs,f_obs, color=:black, ms=2, legend=:none)
  plot!(lambda_obs,f_filtered,color=:black)
  plot!(lambda_obs,continuum,color=:cyan)
  scatter!(lambda_obs[anchors], f_filtered[anchors], color=:blue, ms=4)


if write_output
  jldopen(continuum_ras_fn, "w") do file
    file["order_idx_with_continuum"] = orders_to_use
    file["continuum_anchors_idx"] = anchors
    file["continuum"] = continuum_output
    file["flux_filtered"] = f_filtered
  end
end





using Plots
order_idx = 50
 println("# Order Index: ", order_idx)
 lambda_obs = view(spec.λ,:,order_idx)
 f_obs = view(f_norm,:,order_idx)
 var_obs = view(var_norm,:,order_idx)
 (anchors, continuum, f_filtered_internal) = Continuum.calc_continuum(lambda_obs,f_obs,var_obs,
      stretch_factor=10.0, ν=1.0, merging_threshold=0.25, smoothing_half_width = 40, min_R_factor = 200.0, verbose=false)
  f_filtered = f_filtered_internal
  println(" size(anchors) = ", size(anchors), " size(continuum) = ", size(continuum), " size(f_filtered) = ", size(f_filtered) )

  #plt = scatter(lambda_obs,f_obs, color=:black, ms=2, legend=:none)
  plot!(plt,lambda_obs,f_filtered,color=:black)
  plot!(plt,lambda_obs,continuum,color=:cyan)
  #=
  if length(anchors) >=1
    scatter!(plt,lambda_obs[anchors], f_filtered[anchors], color=:blue, ms=4)
  end
  plt
  =#
  #=
  (anchors2, continuum2, f_filtered_internal2) = Continuum.calc_continuum(lambda,f_norm[:,order_idx],
       stretch_factor=5.0, ν=1.0, merging_threshold=0.25, smoothing_half_width = 40, verbose=false)
  println("# num_anchors1 = ", length(anchors), "  num_anchors2 = ", length(anchors2))
  plot!(lambda,continuum2,color=:magenta)
  scatter!(lambda[anchors2], f_filtered_internal2[anchors2], color=:red, ms=4)
  =#

xmin = first(lambda)-0.05; xmax = xmin+20;
 xlims!(xmin,xmax)
 idx_pix = searchsortedfirst(lambda,xmin):searchsortedlast(lambda,xmax)
 ylims!(NaNMath.extrema(f_obs[idx_pix]))


xmax = last(lambda)+0.05; xmin = xmax-20;
  xlims!(xmin,xmax)
  idx_pix = searchsortedfirst(lambda,xmin):searchsortedlast(lambda,xmax)
  ylims!(NaNMath.extrema(f_obs[idx_pix]))

  #xmin = 4570; xmax = 4580;
  xmin = 4945; xmax = 4965;
    xlims!(xmin,xmax)
    idx_pix = searchsortedfirst(lambda,xmin):searchsortedlast(lambda,xmax)
    ylims!(NaNMath.extrema(f_obs[idx_pix]))

xlims!(10840,10860)
xlims!(5880,5915)
  #xlims!(xmin,xmax)
xlims!(5350,5385)
xlims!(5385,5400)



xlims!(4950,4960)
true

  xlims!(5875,5905)
  #ylims!(0.95,1.15)
  #xmin = minimum(lambda_plt)
  #xmax = minimum(lambda_plt)
  #xlims!(xmin,xmax)
  #ylims!(extrema(f_filtered[xmin.<=lambda.<=xmax]))



  lambda = spec.λ[:,order_idx]
  f_obs = spec.flux[:,order_idx]
  findall(isnan.(f_obs))

 (anchors, continuum, f_filtered_internal) = Continuum.calc_continuum(lambda,f_obs,
       stretch_factor=5.0, ν=1.0, merging_threshold=0.25, verbose=true)


@time continuum = Continuum.calc_continuum(lambda,f_obs,stretch_factor=5.0, ν=1.0, merging_threshold=0.25, verbose=false)


#using Plots
continuum = zeros(size(spec.flux))
  for order_idx in 4:118
    #order_idx =30
    #possible_pix = get_pixel_range(get_inst(spec),order_idx)
    possible_pix = min_pixel_in_order(get_inst(spec)):max_pixel_in_order(get_inst(spec))
    bad_pix = bad_col_ranges(get_inst(spec),order_idx)
    pix_rng = EchelleInstruments.calc_complement_index_ranges(possible_pix,bad_pix)
    pix = mapreduce(p->collect(p),vcat,pix_rng)
    lambda = spec.λ[pix,order_idx]
    f_obs = spec.flux[pix,order_idx]
    (anchors, continuum_order, f_filtered) = Continuum.calc_continuum(lambda,f_obs,λout=spec.λ[possible_pix,order_idx],
         stretch_factor=10.0, merging_threshold=1.0, verbose=true)
    continuum[possible_pix,order_idx] .= continuum_order
    lambda_plt = spec.λ[possible_pix,order_idx]
    plot(lambda_plt,continuum_order)
    scatter!(lambda[anchors], f_filtered[anchors], color=:blue, ms=4)
  end

plot(spec.λ[:,60],continuum[:,60])
true

plot(lambda, spec.flux[:,order_idx]./continuum, legend=:none)
 xlabel!("λ")
 ylabel!("Flux (continuum-normalized)")
 ylims!(0.9,1.05)
 #ylims!(extrema( f_obs./continuum))
 #xlims!(minimum(lambda),5340)
true
=#

#=
  #plot!(lambda,f_clean,ms=1.5, color=:green)
  #plot!(lambda,f_clip_threshold,color=:blue)
  #ylims!(0.9,1.2)
  #plot!(lambda,f_median_filtered.-f_clean)

 #plot!(lambda,rollingpin_radius./rollingpin_R_max)
 #plot(lambda,f_obs, legend=:none)
 #plot(lambda,f_obs, legend=:none)
 #plot!(lambda,f_interp)
 anchor_vals = f_median_filtered[anch]
 #scatter!(lambda[anch],anchor_vals, color=:blue)
 scatter!(lambda[anch],anchor_vals, color=:red)
 scatter!(lambda[anch[anch_mask]],anchor_vals[anch_mask], color=:blue)
 #interp_cubic = CubicSpline(f_median_filtered[anch],anchor_vals )
 #plot!(lambda,interp_cubic.(lambda), color=:green)
 ylims!(0.9,1.2)
 #xlims!(5006.5,5008)

  scatter!(lambda[anchors_merged],f_median_filtered[anchors_merged], color=:magenta)
 plot!(lambda,calc_continuum_from_anchors(lambda,f_median_filtered,anchors_merged),color=:red)

 #=
 plot!(λ_med_filtered,f_median_filtered)
 plot!(λ_med_filtered,f_median_filtered.+iqr)
 plot!(λ_med_filtered,f_median_filtered.-iqr)
 plot!(λ_med_filtered,f_interp, ms=2.0)
 scatter!(λ_med_filtered[anch],f_median_filtered[anch], ms=4.0)
 =#
 #ylims!(0.9,1.03)

#interp_bsp = BSplineApprox(lambda, f_obs,3,4,:ArcLen,:Average)
=#
