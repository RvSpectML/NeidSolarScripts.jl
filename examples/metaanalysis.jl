if occursin(r"RvSpectMLEcoSystem$", pwd())
 cd("NeidSolarScripts")
 using Pkg
 Pkg.activate(".")
end

using CSV, DataFrames, Query
 using JLD2, FileIO
 using Statistics, Dates
 using RvSpectMLBase, RvSpectML
 using Plots

data_to_read = ["times", "rvs_ccf_espresso", "σ_rvs_espresso",
                "times_binned", "rvs_binned", "rv_rms",
                "rv_0", "rv_mean", "rv_slope",
                "order_snr", "order_rvs_g",
                "v_grid", "ccfs_espresso", "ccf_vars_espresso"]
 df = DataFrame(:filename=>readdir("output")) |> @filter(endswith(_.filename,".jld2")) |> DataFrame

 for field in data_to_read
   df[!,Symbol(field)] = map(i->load(joinpath("output",i.filename),field),eachrow(df))
 end
 df[!,"jd"] = map(i->datetime2julian(DateTime(parse(Int,SubString(df.filename[i],7:10)), parse(Int,SubString(df.filename[i],11:12)),
               parse(Int,SubString(df.filename[i],13:14)))), 1:size(df,1) )
 #df = df |> @filter( ! (2459213.4 < _.jd < 2459213.6) ) |> DataFrame
 #df = df |> @filter( ! (2459224.4 < _.jd < 2459224.6) ) |> DataFrame
 t0 = datetime2julian(DateTime(2021,1,22))

plt = scatter(df[!,"jd"].-t0, first.(df[!,"rv_0"]), xlabel="Day", ylabel="RV fit at noon (m/s)", label=:none)
 savefig("rv_solar_noon_vs_t_all.png")

df = df |> @filter( ! (2459213.4 < _.jd < 2459213.6) ) |> DataFrame
  df = df |> @filter( ! (2459224.4 < _.jd < 2459224.6) ) |> DataFrame
  plt = scatter(df[!,"jd"].-t0, first.(df[!,"rv_0"]), xlabel="Day", ylabel="RV fit at noon (m/s)", label=:none)
  first_obs_to_fit = 3
  x = [df[!,"jd"][first_obs_to_fit:end].-t0  ones(length(df[!,"jd"][first_obs_to_fit:end])) ]
  y = first.(df[!,"rv_0"])[first_obs_to_fit:end]
  linfit = (x'*x)\x'*y
  y_fit = (x*linfit)
  slope = linfit[1]
  rms = std(y.-y_fit)
  plot!(plt,x[:,1],y_fit, label=:none)
  title!(plt,"slope = " * string(round(slope,digits=2)) * " m/s/day  rms = " * string(round(rms,digits=2)) * " m/s")
  savefig("rv_solar_noon_vs_t.png")
  #ylims!(plt, 132,143)



plt = scatter(df[!,"jd"].-t0, first.(df[!,"rv_slope"]), xlabel="Day", ylabel="RV slope (m/s/hr)", label=:none)
 savefig("rv_slope_vs_t.png")

plt1 = plot()
 plt2 = plot()
 for i in 1:size(df,1)
    plot!(plt1,df[i,"times"].*24, df[i,"order_snr"][60,:], label=string(df[i,"jd"]-t0), color=i, legend=:none)
    plot!(plt2,df[i,"times"].*24, df[i,"order_snr"][20,:]./df[i,"order_snr"][80,:] .- maximum(df[i,"order_snr"][20,:]./df[i,"order_snr"][80,:]), color=i, label=string(df[i,"jd"]-t0), legend=:none)
 end
 ylims!(plt1,2.8e4,3.2e4)
 plt = plot(plt1,plt2, layout = (2,1))
 xlabel!(plt,"Time (hr - solar noon)")
 ylabel!(plt1,"SNR in Order 60")
 display(plt1)
 savefig("snr_vs_t_day.png")
 ylabel!(plt2,"Δ (SNR_20/SNR_80)")
 display(plt)
 savefig("snr_color_vs_t_day.png")


daily_ccfs = zeros(size(df[1,"ccfs_espresso"],1),size(df,1))
for i in 1:size(daily_ccfs,2)
  idx_to_use = -1.5 .<= df[i,"times"]*24.0 .<= 1.5
  daily_ccfs[:,i] .= vec(mean(df[i,"ccfs_espresso"][:,idx_to_use],dims=2))
end


obsid_start = 3
obsid_stop = size(daily_ccfs,2)-1
 daily_ccfs_norm = daily_ccfs[:,obsid_start:obsid_stop]./maximum(daily_ccfs[:,obsid_start:obsid_stop],dims=1)
 ccf_template = vec(mean(daily_ccfs_norm,dims=2))
 plot(daily_ccfs_norm.-ccf_template)
 heatmap(df[1,"v_grid"],1:length(df[!,"jd"][obsid_start:obsid_stop]),daily_ccfs_norm'.-ccf_template')
 title!("Residuals of Mid-day mean CCF")
 xlabel!("v (m/s)")
 ylabel!("Obs ID")
 xlims!(-30e3,30e3)
 savefig("daily_ccf_resid.png")
