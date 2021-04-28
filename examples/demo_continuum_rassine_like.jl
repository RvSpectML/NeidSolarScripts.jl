using NeidSolarScripts

using Statistics
# Test model/dataset
f_line(x;λo::Real, σ::Real, d::Real) = 1-d*exp(-0.5*((x-λo)/σ)^2)
 function f_model(x; fwhm::Real)
  speed_of_light_mks = 3e8 # m/s
  λcs = [5001, 5002.5, 5003,5004,5008,5008.3]
  ds = [0.25, 0.75, 0.5, 0.5, 0.5, 0.5]
  σs = λcs.*(fwhm/speed_of_light_mks)/sqrt(8(log(2)))
  f = prod(map(i->f_line(x,λo=λcs[i],σ=σs[i],d=ds[i]),1:length(λcs)))
  f *= 1.1-0.1*((x-mean(λcs))/(maximum(λcs)-minimum(λcs)))^2
end

λlo = 5000.0
 λhi = λlo*(5010/5000)
 fwhm_sol = 7.3e3 # m/s
 speed_of_light_mks = 3e8 # m/s
 R = 120000
 snr_per_res = 300
 oversample = 4
 n_res_elements = log(λhi/λlo)*R
 n_pixels = floor(Int64,n_res_elements*oversample/2)*2
 snr_per_pixel = snr_per_res/sqrt(n_pixels/n_res_elements)

 lambda = collect(range(λlo,stop=λhi,length=n_pixels))
 f_true = f_model.(lambda; fwhm=fwhm_sol)
 f_obs = f_true .+ (1/snr_per_pixel).*randn(length(f_true))

#using
(anchors, continuum, f_filtered) = Continuum.calc_continuum(lambda,f_obs)

using Plots
scatter(lambda,f_obs, color=:black, ms=2, legend=:none)
  plot!(lambda,f_filtered,color=:cyan)
  plot!(lambda,continuum,color=:magenta)
  scatter!(lambda[anchors], f_filtered[anchors], color=:blue, ms=4)
  ylims!(0.95,1.15)

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
