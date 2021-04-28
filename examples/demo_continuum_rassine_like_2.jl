using NeidSolarScripts
using EchelleInstruments

fits_path = "/mnt/data_simons/NEID/DRPv0.7-fixedflatfielding2/"
fits_fn = "neidL1_20210104T182500.fits"
continuum_afs_path = joinpath(fits_path,"output","continuum")
continuum_afs_fn = "neidL1_20210104T182500_continuum=afs.jld2"

spec = NEID.read_data(joinpath(fits_path,fits_fn))

#using Plots

order_idx = 60
 lambda = spec.λ[:,order_idx]
 f_obs = spec.flux[:,order_idx]
 (anchors, continuum, f_filtered_internal) = Continuum.calc_continuum(lambda,f_obs,
      stretch_factor=0.5, merging_threshold=0.04, verbose=true)
  println("# num_anchors = ", length(anchors))
 f_filtered = f_filtered_internal
 #f_filtered = Continuum.calc_rolling_median(lambda,f_obs, width=4)
 scatter(lambda,f_obs, color=:black, ms=2, legend=:none)
  plot!(lambda,f_filtered,color=:cyan)
  plot!(lambda,continuum,color=:magenta)
  scatter!(lambda[anchors], f_filtered[anchors], color=:blue, ms=3)
  #ylims!(0.95,1.15)
  xmin = 5360
  xmax = xmin+15
  xlims!(xmin,xmax)
  ylims!(extrema(f_filtered[xmin.<=lambda.<=xmax]))

plot(lambda, f_obs./continuum, legend=:none)
 xlabel!("λ")
 ylabel!("Flux (continuum-normalized)")
 ylims!(0.9,1.05)
 #ylims!(extrema( f_obs./continuum))
 #xlims!(minimum(lambda),5340)
true
true

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
