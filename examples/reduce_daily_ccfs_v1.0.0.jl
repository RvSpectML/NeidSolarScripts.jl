### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ cd026a7a-39f7-429e-95ac-a3be255f1f29
begin
	import Pkg
	Pkg.status()
end

# ╔═╡ 37521853-b548-44b1-b4ea-df73eb2682d2
begin
	using NeidSolarScripts
	using RvSpectMLBase
	using EchelleInstruments
	using EchelleCCFs
	using RvSpectML
	using DataFrames, Query
	using CSV, JLD2, FileIO
	using Dates
	using PlutoUI
	using Plots, ColorSchemes
	using Statistics, NaNMath
	using Polynomials
	using LinearAlgebra, BandedMatrices
	using QuadGK
end

# ╔═╡ 87597e02-19c7-418a-8eb3-ad25695bb5bf
md"""
# Process Daily NEID Solar CCFs
"""

# ╔═╡ ba3a32ee-9664-4279-b27b-0fed992e4d1c
md"## Code "

# ╔═╡ f64f9468-25ee-4799-b776-b1cfbe5b7f57
md"""
### Parameters for processing
"""

# ╔═╡ c259acbf-da3d-4b73-af2b-8fc22f2548d7
begin
	ccf_path = "/mnt/data_simons/NEID/outputs_1.0.0/"
	daily_filename = "daily_ccfs_1.jld2"
	max_days_to_use = 256
	harspn_path = "../data/"
	md"Set which directory/files to process."
end

# ╔═╡ 15741d2d-5fc2-4dc0-baf1-3ba3a72b5399
begin
	solar_hour_angle_threshold_good = 1.5  #2.0
	pyrhelio_ratio_min_good = 0.90
	pyrhelio_ratio_max_good = 1.05
	pyrhelio_rms_min_good = 0.0
	pyrhelio_rms_max_good = 0.003
	min_obs_in_day = 20
	min_obs_binned_in_day = 10 # 20
	md"Set parameters for which observations/days to use for calculations."
end;

# ╔═╡ d6219a57-7c1d-4b23-9b80-ae96f005b762
begin   # Parameters for building CCF template
	solar_hour_angle_threshold_template = 1.0
	pyrhelio_ratio_min_template = 0.95
	pyrhelio_ratio_max_template = 1.05
	pyrhelio_rms_min_template = 0.0000001
	pyrhelio_rms_max_template = 0.002
	md"Set parameters for building template CCF"
end

# ╔═╡ dfea1ed2-4431-4d5b-99c7-560dc4472ab1
begin 
	frac_of_width_to_fit = 1.0  # 1.0
	measure_width_at_frac_depth = 0.5 
	md"Set parameters for measuring RV from CCF"
end

# ╔═╡ 011c5463-2edf-4903-b512-0df12d8b3779
save_figs = false

# ╔═╡ 7f12ff80-e32c-400b-8cdb-1a800fdbc791
md"## Results Summary"

# ╔═╡ 60559545-8f85-4027-8471-4ce3d7cc952a
md"""
### Script
"""

# ╔═╡ d2543921-868a-47d5-94f8-a965161ef8b4
md"### Compute RVs"

# ╔═╡ 474d1a20-91ec-4b50-9caf-1991de12fb9e
function calc_ccf_order_depth(vel::AV1, flux::AV2, var::AV3; smooth_factor::Real = 2) where {T1<:Real, T2<:Real, T3<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3}}
	(f_mean, f_deriv) = RvSpectML.TemporalGPInterpolation.predict_mean_and_deriv(vel, flux, vel;sigmasq_obs=var, use_logx=false, use_logy=false, smooth_factor=smooth_factor)
	depth = 1.0 .- minimum(f_mean)/maximum(f_mean)
end

# ╔═╡ 8d98f83f-cdb6-4856-b9a2-0e8afe39702f
function calc_ccf_order_σ_rv(vel::AV1, flux::AV2, var::AV3; smooth_factor::Real = 2) where {T1<:Real, T2<:Real, T3<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3}}
	(f_mean, f_deriv) = RvSpectML.TemporalGPInterpolation.predict_mean_and_deriv(vel, flux, vel;sigmasq_obs=var, use_logx=false, use_logy=false, smooth_factor=smooth_factor)
	
   exp_sigma_rv = 1 / sqrt( sum(f_deriv.^2 ./ var) )
end

# ╔═╡ a25e4ec0-e8f2-436e-9f3a-c0c7da6f7fe2
#scatter(1:num_order_idx,map(ordix->calc_ccf_order_depth(x["v_grid"],x["order_ccfs"][:,ordix,100], x["order_ccf_vars"][:,6,100]),1:num_order_idx))

# ╔═╡ b03aeea2-4591-4d2a-af22-441e9bb2b571
#= begin
	scatter(1:num_order_idx,map(ordix->calc_ccf_order_σ_rv(x["v_grid"],x["order_ccfs"][:,ordix,100], x["order_ccf_vars"][:,6,100]),1:num_order_idx))
	#ylims!(0,1e6)
end =#

# ╔═╡ 8acc524d-692e-47f4-8c4d-ecbe6a8f4226
function calc_ccf_order_σ_rv_expr(vel::AV1, flux::AV2, var::AV3; smooth_factor::Real = 2) where {T1<:Real, T2<:Real, T3<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3}}
	(f_mean, f_deriv) = RvSpectML.TemporalGPInterpolation.predict_mean_and_deriv(vel, flux, vel;sigmasq_obs=var, use_logx=false, use_logy=false, smooth_factor=smooth_factor)
	lsf_v_width = 3000.0 * sqrt(2)
	pixel_v_width = 620.953 
	mask_v_width = pixel_v_width#lsf_v_width
   covar = zeros(length(f_deriv),length(f_deriv))
	for i in 1:length(vel)
		covar[i,i] = var[i]*1.01
		
		for j in (i+1):length(vel)
			dv = abs(vel[i]-vel[j])
			
			if dv < pixel_v_width
				covar[i,j] = (1-dv/mask_v_width)/mask_v_width*sqrt(var[i]*var[j])
			end
			covar[i,j] += sqrt(var[i]*var[j])*exp(-0.5*(dv/lsf_v_width)^2)/(lsf_v_width*sqrt(2π))  
			covar[j,i] = covar[i,j]
		end	
	end
   exp_sigma_rv = 1 / sqrt( sum( (f_deriv' * (covar \ f_deriv) ) ) ) 
end

# ╔═╡ 493fea80-ad5a-4fb1-add5-ad81e262ea9f
function calc_ccf_order_σ_rv_expr2(vel::AV1, flux::AV2, var::AV3; smooth_factor::Real = 2) where {T1<:Real, T2<:Real, T3<:Real, AV1<:AbstractVector{T1}, AV2<:AbstractVector{T2}, AV3<:AbstractVector{T3}}
	(f_mean, f_deriv) = RvSpectML.TemporalGPInterpolation.predict_mean_and_deriv(vel, flux, vel;sigmasq_obs=var, use_logx=false, use_logy=false, smooth_factor=smooth_factor)
	num_vels = length(vel)
	pixel_v_width = 620.953 
	lsf_v_width = 3000.0 
	lsf_σ = lsf_v_width /sqrt(8*log(2))  
	mask_v_width =  lsf_v_width #/sqrt(8*log(2)) # lsf_v_width # 
	δv = vel[2]-vel[1]
	Δv_off_dag = 3*max(lsf_v_width,mask_v_width)
	ndiagonals = searchsortedfirst(vel,first(vel)+Δv_off_dag) -1
	if ndiagonals < length(vel)
		covar = BandedMatrix{Float64}(undef, (num_vels,num_vels), (ndiagonals,ndiagonals))
		covar[band(0)] .= 1.00000001 # var
		denom = quadgk(x->max((1-0/mask_v_width)/mask_v_width,0.0)*exp(-0.5*((x-0)/lsf_σ)^2)/(lsf_σ*sqrt(2π)),-6*lsf_σ,6*lsf_σ)[1]*δv
		#denom = quadgk(x->exp(-0.5*(0/mask_v_width)^2)/(mask_v_width*sqrt(2π))*exp(-0.5*((x-0)/lsf_σ)^2)/(lsf_σ*sqrt(2π)),-6*lsf_σ,+6*lsf_σ)[1]*δv
		#denom = (1/(lsf_σ*sqrt(2π)))
		for i in 1:ndiagonals
			dv = vel[1+i]-vel[1]
			#val = max(1-dv/pixel_v_width,0.0) 
			#val += exp(-0.5*(dv/lsf_σ)^2)  
			val = quadgk(x->max((1-dv/mask_v_width)/mask_v_width,0.0)*exp(-0.5*((x-dv)/lsf_σ)^2/2)/(lsf_σ*sqrt(2π*2)),dv-6*lsf_σ,dv+6*lsf_σ)[1]*δv/denom
			#val = quadgk(x->exp(-0.5*(dv/mask_v_width)^2)/(mask_v_width*sqrt(2π))*exp(-0.5*((x-dv)/lsf_σ)^2)/(lsf_σ*sqrt(2π)),dv-6*lsf_σ,dv+6*lsf_σ)[1]*δv/denom
			#val = (exp(-0.5*(dv/lsf_σ)^2)/(lsf_σ*sqrt(2π)))/denom
			covar[band(i)] .= covar[band(-i)] .= val
		end
		for i in 1:num_vels
			covar[i,BandRange] .*= sqrt.(var[i])
			covar[BandRange,i] .*= sqrt.(var[i])
		end		
		#covar = SymmetricBanded(covar)
	else 
		covar = Symmetric(diagm(var))
	end
   exp_sigma_rv = 1 / sqrt( sum( (f_deriv' * (covar \ f_deriv) ) ) ) 
end

# ╔═╡ 7c302268-c603-4e46-adb1-ba047429c2e0
#calc_ccf_order_σ_rv_expr2(x["v_grid"],x["order_ccfs"][:,6,100], x["order_ccf_vars"][:,6,100],smooth_factor=2)
3000 /sqrt(8*log(2)) 

# ╔═╡ 3f1b1f9c-2ee0-4c8b-832f-d984669450ba
#= begin 
	orders_plt = 1:14 # num_order_idx
	local depth = map(ordix->calc_ccf_order_depth(x["v_grid"],x["order_ccfs"][:,ordix,100], x["order_ccf_vars"][:,ordix,100]),orders_plt)
	scatter(orders_plt,depth)
	local snr = map(ordix->mean(x["order_ccfs"][:,ordix,100]./sqrt.(x["order_ccf_vars"][:,ordix,100]),dims=1)[1],orders_plt)
	#scatter(orders_plt,snr)
	local sigmarv = map(ordix->calc_ccf_order_σ_rv(x["v_grid"],x["order_ccfs"][:,ordix,100], x["order_ccf_vars"][:,ordix,100]),orders_plt)
	#scatter(depth.*sqrt.(snr),sigmarv, label="Diagonal terms only")
	scatter(orders_plt, sigmarv, label="Diagonal terms only")
	local sigmarv2 = map(ordix->calc_ccf_order_σ_rv_expr(x["v_grid"],x["order_ccfs"][:,ordix,100], x["order_ccf_vars"][:,ordix,100]),orders_plt)
	#scatter!(depth.*sqrt.(snr),sigmarv2, ms=2.5,label="Including Off-diagonal terms")
	local sigmarv3 = map(ordix->calc_ccf_order_σ_rv_expr2(x["v_grid"],x["order_ccfs"][:,ordix,100], x["order_ccf_vars"][:,ordix,100],smooth_factor=2),orders_plt)
	#scatter!(depth.*sqrt.(snr),sigmarv3, ms=2.5,label="Experimental")
	scatter!(orders_plt,sigmarv3, ms=2.5,label="Experimental")
	scatter!(orders_plt,day_orders.rms_order_rv_binned[71])	
	#xlims!(20,maximum(depth.*sqrt.(snr)))
	#ylims!(0,10)
	#scatter!(depth.*sqrt.(snr),(sigmarv2./sigmarv), ms=2.5,label="Ratio Sq 2")
	#scatter!(depth.*sqrt.(snr),(sigmarv3./sigmarv), ms=2.5,label="Ratio Sq 3")
	ylims!(0,2)
end =#

# ╔═╡ 58086727-ef4a-4456-bfdf-f820432f4cad
# x["order_ccfs"][300,6,100], x["order_ccf_vars"][300,6,100]

# ╔═╡ 2a41bd71-3796-4a21-b7ce-da0b09d043de
begin
	harpsn_7day = CSV.read(joinpath("..","data","rms7.csv"),DataFrame) 
	rename!(harpsn_7day, [Symbol("# JD bin center"),Symbol("rv rms in bin"), Symbol("number of spec in bin")] .=>  [:jd, :rms_rv, :num_days])
	harpsn_7day = harpsn_7day |> @filter( _.num_days >= 5) |> DataFrame

	harpsn_30day = CSV.read(joinpath("..","data","rms30.csv"),DataFrame) 
	rename!(harpsn_30day, [Symbol("# JD bin center"),Symbol("rv rms in bin"), Symbol("number of spec in bin")] .=>  [:jd, :rms_rv, :num_days])
	harpsn_30day = harpsn_30day |> @filter( _.num_days >= 10) |> DataFrame
end;

# ╔═╡ 21d22e3c-a56a-4c03-a9aa-226acea40c72
md"""
# Results
## One Day
"""

# ╔═╡ fe0b7b47-39bb-45f8-873e-9a4a07d7d1b3
md"""
Order RV slope limits (lo, hi) $(@bind order_slope_lo Slider(-100.0:0.0; default=-40))  $(@bind order_slope_hi Slider(0.0:100.0; default=40.0))

"""

# ╔═╡ b6da427f-0c72-4de7-bdca-b85688ac9bed
md"""
## Combined RVs
#### Order indices to include in RV calculation
1 $(@bind use_ord_1 CheckBox(default=true))
2 $(@bind use_ord_2 CheckBox(default=true))
3 $(@bind use_ord_3 CheckBox(default=true))
4 $(@bind use_ord_4 CheckBox(default=true))
5 $(@bind use_ord_5 CheckBox(default=true))
6 $(@bind use_ord_6 CheckBox(default=true))
7 $(@bind use_ord_7 CheckBox(default=false))
8 $(@bind use_ord_8 CheckBox(default=true))
9 $(@bind use_ord_9 CheckBox(default=true))
10 $(@bind use_ord_10 CheckBox(default=true))
11 $(@bind use_ord_11 CheckBox(default=false))
12 $(@bind use_ord_12 CheckBox(default=true))
13 $(@bind use_ord_13 CheckBox(default=true))
14 $(@bind use_ord_14 CheckBox(default=true))
15 $(@bind use_ord_15 CheckBox(default=false))
16 $(@bind use_ord_16 CheckBox(default=true))
17 $(@bind use_ord_17 CheckBox(default=true))
18 $(@bind use_ord_18 CheckBox(default=false))
19 $(@bind use_ord_19 CheckBox(default=false))
20 $(@bind use_ord_20 CheckBox(default=true))

21 $(@bind use_ord_21 CheckBox(default=true))
22 $(@bind use_ord_22 CheckBox(default=false))
23 $(@bind use_ord_23 CheckBox(default=false))
24 $(@bind use_ord_24 CheckBox(default=false))
25 $(@bind use_ord_25 CheckBox(default=false))
26 $(@bind use_ord_26 CheckBox(default=false))
27 $(@bind use_ord_27 CheckBox(default=false))
28 $(@bind use_ord_28 CheckBox(default=false))
29 $(@bind use_ord_29 CheckBox(default=false))
30 $(@bind use_ord_30 CheckBox(default=false))
31 $(@bind use_ord_31 CheckBox(default=false))
32 $(@bind use_ord_32 CheckBox(default=false))
33 $(@bind use_ord_33 CheckBox(default=false))
34 $(@bind use_ord_34 CheckBox(default=false))
35 $(@bind use_ord_35 CheckBox(default=false))
36 $(@bind use_ord_36 CheckBox(default=false))
37 $(@bind use_ord_37 CheckBox(default=false))
38 $(@bind use_ord_38 CheckBox(default=false))
39 $(@bind use_ord_39 CheckBox(default=false))
40 $(@bind use_ord_40 CheckBox(default=false))

41 $(@bind use_ord_41 CheckBox(default=false))
42 $(@bind use_ord_42 CheckBox(default=false))
43 $(@bind use_ord_43 CheckBox(default=false))
44 $(@bind use_ord_44 CheckBox(default=false))
45 $(@bind use_ord_45 CheckBox(default=false))
46 $(@bind use_ord_46 CheckBox(default=false))
47 $(@bind use_ord_47 CheckBox(default=false))
48 $(@bind use_ord_48 CheckBox(default=false))
49 $(@bind use_ord_49 CheckBox(default=false))
50 $(@bind use_ord_50 CheckBox(default=false))
51 $(@bind use_ord_51 CheckBox(default=false))
52 $(@bind use_ord_52 CheckBox(default=false))
53 $(@bind use_ord_53 CheckBox(default=false))
Manual override $(@bind use_ord_manual CheckBox(default=false))
Highlight $(@bind order_to_highlight NumberField(0:53; default=8))  
"""

# ╔═╡ b94f222a-0c35-49d7-94f2-a2347efb6b2d
#order_idx_for_rvs_manual = vcat(1,4:6,10:11,19:23);
#order_idx_for_rvs_manual = vcat(1,4:6,8,10:12,15:17,19:21);
#order_idx_for_rvs_manual = vcat(13:14,16:21,23,24:25);
#order_idx_for_rvs_manual = vcat(37:41,46:50);
#order_idx_for_rvs_manual = vcat(1:6,8:10,12);
order_idx_for_rvs_manual = vcat(1:14,17:21) #,23:29)
#order_idx_for_rvs_manual = vcat(18:21,23:29)
#order_idx_for_rvs_manual = vcat(36:43,45:50)

# ╔═╡ e7541e93-8b7a-4c26-901e-d887931b8103
begin
	if use_ord_manual
		order_idx_for_rvs = order_idx_for_rvs_manual
	else
		order_idx_for_rvs = Int64[] #order_idx_for_rvs
	if use_ord_1 push!(order_idx_for_rvs,1) end
	if use_ord_2 push!(order_idx_for_rvs,2) end
	if use_ord_3 push!(order_idx_for_rvs,3) end
	if use_ord_4 push!(order_idx_for_rvs,4) end
	if use_ord_5 push!(order_idx_for_rvs,5) end
	if use_ord_6 push!(order_idx_for_rvs,6) end
	if use_ord_7 push!(order_idx_for_rvs,7) end
	if use_ord_8 push!(order_idx_for_rvs,8) end
	if use_ord_9 push!(order_idx_for_rvs,9) end
	if use_ord_10 push!(order_idx_for_rvs,10) end
	if use_ord_11 push!(order_idx_for_rvs,11) end
	if use_ord_12 push!(order_idx_for_rvs,12) end
	if use_ord_13 push!(order_idx_for_rvs,13) end
	if use_ord_14 push!(order_idx_for_rvs,14) end
	if use_ord_15 push!(order_idx_for_rvs,15) end
	if use_ord_16 push!(order_idx_for_rvs,16) end
	if use_ord_17 push!(order_idx_for_rvs,17) end
	if use_ord_18 push!(order_idx_for_rvs,18) end
	if use_ord_19 push!(order_idx_for_rvs,19) end
	if use_ord_20 push!(order_idx_for_rvs,20) end
	if use_ord_21 push!(order_idx_for_rvs,21) end
	if use_ord_22 push!(order_idx_for_rvs,22) end
	if use_ord_23 push!(order_idx_for_rvs,23) end
	if use_ord_24 push!(order_idx_for_rvs,24) end
	if use_ord_25 push!(order_idx_for_rvs,25) end
	if use_ord_26 push!(order_idx_for_rvs,26) end
	if use_ord_27 push!(order_idx_for_rvs,27) end
	if use_ord_28 push!(order_idx_for_rvs,28) end
	if use_ord_29 push!(order_idx_for_rvs,29) end
	if use_ord_30 push!(order_idx_for_rvs,30) end
	if use_ord_31 push!(order_idx_for_rvs,31) end
	if use_ord_32 push!(order_idx_for_rvs,32) end
	if use_ord_33 push!(order_idx_for_rvs,33) end
	if use_ord_34 push!(order_idx_for_rvs,34) end
	if use_ord_35 push!(order_idx_for_rvs,35) end
	if use_ord_36 push!(order_idx_for_rvs,36) end
	if use_ord_37 push!(order_idx_for_rvs,37) end
	if use_ord_38 push!(order_idx_for_rvs,38) end
	if use_ord_39 push!(order_idx_for_rvs,39) end
	if use_ord_40 push!(order_idx_for_rvs,40) end
	if use_ord_41 push!(order_idx_for_rvs,41) end
	if use_ord_42 push!(order_idx_for_rvs,42) end
	if use_ord_43 push!(order_idx_for_rvs,43) end
	if use_ord_44 push!(order_idx_for_rvs,44) end
	if use_ord_45 push!(order_idx_for_rvs,45) end
	if use_ord_46 push!(order_idx_for_rvs,46) end
	if use_ord_47 push!(order_idx_for_rvs,47) end
	if use_ord_48 push!(order_idx_for_rvs,48) end
	if use_ord_49 push!(order_idx_for_rvs,49) end
	if use_ord_50 push!(order_idx_for_rvs,50) end
	if use_ord_51 push!(order_idx_for_rvs,51) end
	if use_ord_52 push!(order_idx_for_rvs,52) end
	if use_ord_53 push!(order_idx_for_rvs,53) end		
	end
	order_idx_for_rvs
end;

# ╔═╡ 473302a3-2541-49a1-8f21-906bc2c9d279
md"""
Order RV RMS limits (lo, hi) 
$(@bind order_rv_rms_lo Slider(0.0:0.1:2.0; default=0))  
$(@bind order_rv_rms_hi Slider(0.5:0.2:200.0; default=6))
"""

# ╔═╡ 9833f763-ece1-4739-a6d1-0f3182d32d5d
md"### Checking for Chromatic Effects"

# ╔═╡ 98087935-894b-408c-b708-02f05d97a96c
#f = JLD2.jldopen(joinpath(ccf_path,"20210208/daily_ccfs_7.jld2"))

# ╔═╡ 2c88f409-6763-41d8-b1cc-168df855728b
md"""
# Implementation Details
"""

# ╔═╡ 38bb2810-3a30-4a39-99dd-83518c093be2
md"""
## Functions used above
"""

# ╔═╡ c23c4063-947e-436a-910b-c43afdb9c244
function make_dataframe_of_ccf_files(path::String,daily_filename::String)
	datafiles = DataFrame(:datestr=>String[],:filename=>String[])
	for (root, dirs, files) in walkdir(path)
	    for file in files
			fn = joinpath(root, file)
			if !(isfile(fn) && filesize(fn)>0) continue end
			m = match(r"(\d+)/" * daily_filename,fn)
			if m == nothing continue end
			datestr = m.captures[1]
			d = Dict(:datestr=>datestr, :filename=>fn)
			push!(datafiles,d)
	    end
	end
	return datafiles
end

# ╔═╡ 3c4b1028-bcf7-45b0-bba2-dab5a63724d0
datafiles_orig = make_dataframe_of_ccf_files(ccf_path,daily_filename);

# ╔═╡ f83d3db6-39bb-4f3c-9870-790b37741646
begin
	df_wavecal_quality = CSV.read(joinpath(ccf_path,"..","Wavecal_Solar_Days_Quality.csv"),DataFrame)
	df_wavecal_quality[!,:datestr] = Dates.format.(df_wavecal_quality[!,Symbol("UT Day")],"YYYYmmdd") 
	#delete!(df_wavecal_quality,Symbol("UT Day"))
	rename!(df_wavecal_quality,	Symbol("LFC flag")=>:lfc_flag)
	rename!(df_wavecal_quality,	Symbol("Drift flag")=>:drift_flag)
	df_wavecal_quality[!,:lfc_flag] = convert.(Bool,df_wavecal_quality[!,:lfc_flag])
	df_wavecal_quality[!,:drift_flag] =convert.(Bool,df_wavecal_quality[!,:drift_flag])
	#df_wavecal_quality 
	datafiles = innerjoin(datafiles_orig,df_wavecal_quality,on = :datestr) 
end;

# ╔═╡ 26caf850-c46e-453c-aad9-c7bdb71f858b
function select_days_to_analyze_add_manifest(df::DataFrame; solar_hour_angle_threshold::Real = solar_hour_angle_threshold_good,  pyrhelio_ratio_min::Real = pyrhelio_ratio_min_good, pyrhelio_ratio_max::Real = pyrhelio_ratio_max_good, pyrhelio_rms_max::Real = pyrhelio_rms_max_good, pyrhelio_rms_min::Real = pyrhelio_rms_min_good, min_obs_in_day::Integer = 5, max_days_to_use::Integer = 1000)
	df[!,:nobs_in_file] .= 0
	df[!,:analyze_day] .= false
	df[!,:manifest] = fill(DataFrame(),size(df,1))
	
	num_days_to_use_found = 0
	for (i,day) in enumerate(eachrow(df))
		if day.lfc_flag || day.drift_flag
			continue
		end
		solar_info_fn = joinpath(dirname(day.filename),"..",day.datestr * "_solar_info.csv")
		println("Looking for ", solar_info_fn)
		if !( isfile(solar_info_fn) && filesize(solar_info_fn)>0)
			println("# Warning missing ", solar_info_fn, " skipping that day.")
			continue
		end
		@assert isfile(solar_info_fn) && filesize(solar_info_fn)>0
		df_solarinfo = CSV.read(solar_info_fn, DataFrame)
		manifest_all = load(day.filename,"manifest")
		df_combo = innerjoin(manifest_all, df_solarinfo, on = :bjd, makeunique = true)
		println("# Sizes: manifest= ", size(manifest_all,1), " solarinfo= ", size(df_solarinfo), " df_combo= ", size(df_combo,1))
		manifest_use = df_combo |> 
			@filter(abs(_.solar_hour_angle)<= solar_hour_angle_threshold) |> 
			@filter( pyrhelio_ratio_min <= _.pyrheliometer_flux_over_model <= pyrhelio_ratio_max ) |>
			@filter( pyrhelio_rms_min <= _.pyrheliometer_rms_flux <= pyrhelio_rms_max ) |>
				DataFrame
		
		day.nobs_in_file = size(manifest_all,1)
		nobs_use = size(manifest_use,1)
		day.manifest = df_combo
		#day.manifest[!,:pyrheliometer_flux_over_model] = manifest_use[!,pyrheliometer_flux_over_model]
		#day.manifest[!,:pyrheliometer_rms_flux] = manifest_use[!,:pyrheliometer_rms_flux]
		if nobs_use <= min_obs_in_day
			continue 
		end
		day.analyze_day = true
		num_days_to_use_found += 1
		if num_days_to_use_found >= max_days_to_use
			break
		end
		
	end
	df_out = df |> @filter(_.analyze_day) |> @orderby_descending(_.datestr)|>  @take(max_days_to_use) |> @orderby(_.datestr) |> DataFrame
 	return df_out
end

# ╔═╡ 785a2042-bb12-4a2a-8117-c2cfb396d29c
datafiles_use = select_days_to_analyze_add_manifest(datafiles,
					solar_hour_angle_threshold=solar_hour_angle_threshold_good,
					pyrhelio_ratio_min=pyrhelio_ratio_min_good, pyrhelio_ratio_max=pyrhelio_ratio_max_good,
					pyrhelio_rms_max=pyrhelio_rms_max_good, pyrhelio_rms_min=pyrhelio_rms_min_good,
					min_obs_in_day=min_obs_in_day, max_days_to_use=max_days_to_use);

# ╔═╡ 0b6ec47b-3610-453f-b0c3-a99d5fb83ab4
begin
	local plt = plot()
	hist1 = Float64[]
	hist2 = Float64[]
	hist3 = Float64[]
	hist4 = Float64[]
	for row in eachrow(datafiles_use)
		if row.drift_flag || row.lfc_flag
			continue
		end
		sha_mask = abs.(row.manifest.solar_hour_angle).<= solar_hour_angle_threshold_good
		pyrh_mask1 = (pyrhelio_ratio_min_good .<= row.manifest.pyrheliometer_flux_over_model .<= pyrhelio_ratio_max_good ) 
		pyrh_mask2 = ( pyrhelio_rms_min_good .< row.manifest.pyrheliometer_rms_flux .<= pyrhelio_rms_max_good )
		mask = sha_mask .& pyrh_mask1  #.& pyrh_mask2
		
		histogram!(plt,row.manifest[mask,:pyrheliometer_rms_flux],bins=range(0,stop=pyrhelio_rms_max_good*2,length=40),alpha=0.3,label=:none)
	end
	xlims!(plt,0,pyrhelio_rms_max_good*2)
	xlabel!(plt,"RMS Pyrheliometer flux")
	ylabel!(plt,"Count")
	title!(plt,"RMS Pyrheliometer flux within selected observations")
	plt
end

# ╔═╡ c1bf8bf8-fb05-4731-8140-c0f6332496bd
begin
	local plt = plot()
	num_passed1 = Int64[]
	num_passed2 = Int64[]
	num_passed3 = Int64[]
	num_passed = Int64[]
	for row in eachrow(datafiles_use)
		 if row.drift_flag || row.lfc_flag
			continue
		end
		sha_mask = abs.(row.manifest.solar_hour_angle).<= solar_hour_angle_threshold_good
		pyrh_mask1 = (pyrhelio_ratio_min_good .<= row.manifest.pyrheliometer_flux_over_model .<= pyrhelio_ratio_max_good ) 
		pyrh_mask2 =  ( pyrhelio_rms_min_good .< row.manifest.pyrheliometer_rms_flux .<= pyrhelio_rms_max_good )
		mask = sha_mask .& pyrh_mask2  .& pyrh_mask1
		push!(num_passed1, sum(sha_mask))
		push!(num_passed2, sum(pyrh_mask1))
		push!(num_passed3, sum(pyrh_mask2))
		push!(num_passed, sum(mask))
		histogram!(plt,row.manifest[mask,:pyrheliometer_flux_over_model],bins=range(max(0.85,pyrhelio_ratio_min_good),stop=min(1.1,pyrhelio_ratio_max_good),length=40),alpha=0.3, label=:none)
	end
	xlims!(plt,max(0.85,pyrhelio_ratio_min_good),min(1.1,pyrhelio_ratio_max_good))
	xlabel!(plt,"Ratio of Pyrheliometer Flux/Model")
	ylabel!(plt,"Count")
	title!(plt,"Distribution of Pyrheliometer Flux/Model")
	plt
end

# ╔═╡ 1ee5e5d4-9b2a-425e-bb3e-386c8db8cd27
begin
	local plt = plot()
	histogram!(plt,num_passed,bins=range(0,maximum(num_passed3),length=40),alpha=0.3,label="Combo")
	histogram!(plt,num_passed1,bins=range(0,maximum(num_passed3),length=40),alpha=0.3,label="SHA")
	histogram!(plt,num_passed2,bins=range(0,maximum(num_passed3),length=40),alpha=0.3,label="Ratio Pyrh")
	histogram!(plt,num_passed3,bins=range(0,maximum(num_passed3),length=40),alpha=0.3,label="RMS Pyrh")
	xlabel!(plt,"Observations per day")
	ylabel!(plt,"Count")
	title!(plt,"Distributions of observations/day passing each cut")
end

# ╔═╡ 0e1f6708-53d8-4cca-acf3-e7b2cac1dd7a
begin   # Get size of order ccfs
	v_grid =  load(first(datafiles_use.filename),"v_grid")
	num_v = length(v_grid)
	ex_order_ccfs = load(first(datafiles_use.filename),"order_ccfs")
	num_order_idx = size(ex_order_ccfs,2)
end;

# ╔═╡ 02e9bc9e-9d86-4b1f-adef-38a67d6ba2d1
v_grid

# ╔═╡ ebf32f98-3ba0-4650-83ca-4d7a5b310c05
md"""
## Orders RVs
Order index (min, max, step) 
$(@bind order_idx_lo NumberField(1:num_order_idx; default=1))  
$(@bind order_idx_hi NumberField(1:num_order_idx; default=12))
$(@bind order_idx_step NumberField(1:10; default=1))
Use order range? $(@bind use_order_range CheckBox(default=true))
Plot bad days? $(@bind plt_all_days CheckBox(default=false))
"""

# ╔═╡ b886c106-b325-464b-8e71-6dc46350e48c
orders_plt = use_order_range ? (order_idx_lo:order_idx_step:order_idx_hi) : order_idx_for_rvs;

# ╔═╡ b957ac39-458f-48b6-8683-229d09ecaea4
function select_days_to_analyze_add_manifest_old(df::DataFrame; solar_hour_angle_threshold::Real = solar_hour_angle_threshold_good, order_snr_threshold::Real = order_snr_threshold_max*0.9, order_idx_to_use_for_snr::Integer = 60, min_obs_in_day::Integer = 5, max_days_to_use::Integer = 10)
	df[!,:nobs_in_file] .= 0
	#df[!,:nobs_use] .= 0
	df[!,:analyze_day] .= false
	df[!,:manifest] = fill(DataFrame(),size(df,1))
	
	num_days_to_use_found = 0
	for (i,day) in enumerate(eachrow(df))
		manifest_all = load(day.filename,"manifest")
		manifest_use = manifest_all |> 
			@filter(abs(_.solar_hour_angle)<= solar_hour_angle_threshold) |> 
				DataFrame
		day.nobs_in_file = size(manifest_all,1)
		#day.
		nobs_use = size(manifest_use,1)
		day.manifest = manifest_all
		if nobs_use <= min_obs_in_day
			continue 
		end
		day.analyze_day = true
		num_days_to_use_found += 1
		if num_days_to_use_found >= max_days_to_use
			break
		end
		
	end
	df_out = df |> @filter(_.analyze_day) |> @orderby_descending(_.datestr)|>  @take(max_days_to_use) |> @orderby(_.datestr) |> DataFrame
 	return df_out
end

# ╔═╡ d9877f46-c461-4615-bb0c-0c7541238899
function compute_mean_order_ccf(df::DataFrame; solar_hour_angle_threshold::Real = solar_hour_angle_template, pyrhelio_ratio_min::Real = pyrhelio_ratio_min_template, pyrhelio_ratio_max::Real = pyrhelio_ratio_max_template, pyrhelio_rms_max::Real = pyrhelio_rms_max_template, pyrhelio_rms_min::Real = pyrhelio_rms_min_template )
	ex_order_ccfs = load(first(df.filename),"order_ccfs")
	mean_order_ccf = zeros(size(ex_order_ccfs,1),size(ex_order_ccfs,2))
	weights_order_ccf = zeros(size(ex_order_ccfs,1),size(ex_order_ccfs,2))
	for row in eachrow(df)
		sha_mask = abs.(row.manifest.solar_hour_angle).<= solar_hour_angle_threshold
		pyrh_mask = (pyrhelio_ratio_min .<= row.manifest.pyrheliometer_flux_over_model .<= pyrhelio_ratio_max ) .& ( pyrhelio_rms_min .< row.manifest.pyrheliometer_rms_flux .<= pyrhelio_rms_max )
	#    snr_mask = map(i->row.manifest.order_snrs[i][order_idx_to_use_for_snr],1:length(row.manifest.order_snrs)) .>= order_snr_threshold
		mask = sha_mask .& pyrh_mask #.& snr_mask
		nobs_use = sum(mask)
		if nobs_use < 60 continue end
		v_grid_this =  load(row.filename,"v_grid")
		@assert all(v_grid .== v_grid_this)
		order_ccfs = load(row.filename,"order_ccfs")
		order_ccf_vars = load(row.filename,"order_ccf_vars")
		@assert size(order_ccfs,1) == size(mean_order_ccf,1)
		@assert size(order_ccfs,2) == size(mean_order_ccf,2)
		num_obs = size(order_ccfs,3)
		weights = ones(num_obs)
		#for j in 1:num_order_idx
		for (j,ord_idx) in enumerate(orders_to_use_default(NEID2D())[1:3])
			if any(order_ccf_vars[:,j,:] .==0)
				println("# Order CCF for ", row.filename, " order idx ", ord_idx, " contains a zero.")
			end
			if any(isnan.(order_ccf_vars[:,j,:]))
				println("# Order CCF for ", row.filename, " order idx ", ord_idx, " contains a NaN.")
			end
			if any(isinf.(order_ccf_vars[:,j,:]))
				println("# Order CCF for ", row.filename, " order idx ", ord_idx, " contains a Inf.")
			end
		end
		for j in 1:num_order_idx
			if any(order_ccf_vars[:,j,mask] .==0) continue end
			#if any(isnan.(order_ccf_vars[:,j,mask]) ) continue end
			for i in 1:num_v
				#weights .= 1.0 ./ order_ccf_vars[i,j,:]
				mean_order_ccf[i,j] += sum(order_ccfs[i,j,mask] .* weights[mask])
				weights_order_ccf[i,j] += sum(weights[mask])
			end
		end
	end
	mean_order_ccf ./= weights_order_ccf
	return mean_order_ccf
end

# ╔═╡ de7e3ebb-f9d2-4830-9bec-f034d965ab76
mean_order_ccf = compute_mean_order_ccf(datafiles_use, solar_hour_angle_threshold=solar_hour_angle_threshold_template, 
	pyrhelio_ratio_min=pyrhelio_ratio_min_template, pyrhelio_ratio_max=pyrhelio_ratio_max_template, pyrhelio_rms_max=pyrhelio_rms_max_template, pyrhelio_rms_min=pyrhelio_rms_min_template )

# ╔═╡ b12592bf-ec54-435f-b795-183cf84bcb4c
begin 
	mrv = []
	for ord_idx in 1:num_order_idx
		templ = vec(mean_order_ccf[:,ord_idx])
		if all(isnan.(templ)) 
			println("# CCF Template for order index ", ord_idx, " is all NaNs.")
			push!(mrv,missing) 
			continue
		end
		
		#deriv = numerical_deriv(v_grid, templ)
		amin = argmin(templ)
		println(ord_idx,": ",amin, " of ", length(v_grid))
		if amin==length(v_grid)
			println(templ)
		end
		mrv_this_order = RVFromCCF.MeasureRvFromCCFTemplate(v_grid=v_grid,template=templ, frac_of_width_to_fit=frac_of_width_to_fit, measure_width_at_frac_depth=measure_width_at_frac_depth)
		push!(mrv,mrv_this_order)
	end
	mrv
end

# ╔═╡ 93465fa0-83bc-413b-b46e-7ed0d3d0d1a7
function calc_all_rvs(fn::String) #, mask::BitVector)
	v_grid_this =  load(fn,"v_grid")
	@assert all(v_grid .== v_grid_this)
	order_ccfs = load(fn,"order_ccfs")
	order_ccf_vars = load(fn,"order_ccf_vars")
	
	num_obs = size(order_ccfs,3) #sum(mask)
	@assert num_obs <= size(order_ccfs,3)
	rv_matrix = zeros(num_order_idx, num_obs)
	σ_rv_matrix = zeros(num_order_idx, num_obs)
	for (i,obs_idx) in enumerate((1:size(order_ccfs,3))) #[mask])
		for ord_idx in 1:num_order_idx
			if typeof(mrv[ord_idx]) == Missing  continue end
			(rv, σ_rv) = mrv[ord_idx](v_grid,order_ccfs[:,ord_idx,obs_idx],order_ccf_vars[:,ord_idx,obs_idx]) 
			rv_matrix[ord_idx,i] = rv
			σ_rv_matrix[ord_idx,i] = σ_rv
		end
	end
	return (rvs=rv_matrix, σ_rvs = σ_rv_matrix, filename=fn)
end

# ╔═╡ 91d634c4-e4ec-49b8-a035-5df0740f9fb0
all_rvs = DataFrame(map(i->calc_all_rvs(datafiles_use.filename[i]),1:length(datafiles_use.filename)));

# ╔═╡ e00fc7cb-05a7-4215-b7f6-5bad647ab090
df_rvs = innerjoin(datafiles_use, all_rvs, on = :filename);

# ╔═╡ 21e0d0a9-0fb4-4bd5-93a3-517481b753d7
(1:size(first(all_rvs.rvs),1))[map(i->any(all_rvs.rvs[1][i,:].==0),1:size(first(all_rvs.rvs),1))]

# ╔═╡ 4226fab5-719e-45a5-ad3d-af712dd743ce
function find_all_good_day_idx(df::DataFrame; solar_hour_angle_threshold::Real = solar_hour_angle_threshold_good, pyrhelio_ratio_min::Real = pyrhelio_ratio_min_good, pyrhelio_ratio_max::Real = pyrhelio_ratio_max_good, pyrhelio_rms_max::Real = pyrhelio_rms_max_good, pyrhelio_rms_min::Real = pyrhelio_rms_min_good )
	local num_good_obs = df |> @map( sum((abs.(_.manifest.solar_hour_angle).<solar_hour_angle_threshold) .&
			(pyrhelio_ratio_min_good .<= _.manifest.pyrheliometer_flux_over_model.<= pyrhelio_ratio_max_good ) .&
			(pyrhelio_rms_min_good .<= _.manifest.pyrheliometer_rms_flux .<= pyrhelio_rms_max_good ) ) 
		) 	|> collect	
	(1:size(df,1))[(num_good_obs .>60) .& .!df.lfc_flag .& .!df.drift_flag 	]
end

# ╔═╡ df6493dd-5730-4993-b45d-495e62f34636
obs_idx_all_good_days = find_all_good_day_idx(datafiles_use)

# ╔═╡ fb55b1e8-eb70-45e9-9139-38ba64a50205
datafiles_use.datestr[obs_idx_all_good_days]

# ╔═╡ 022f7788-a2d0-4489-af0b-71692f800fdf
function find_good_day_idx(df::DataFrame; solar_hour_angle_threshold=solar_hour_angle_threshold_good, pyrhelio_ratio_min::Real = pyrhelio_ratio_min_good, pyrhelio_ratio_max::Real = pyrhelio_ratio_max_good, pyrhelio_rms_max::Real = pyrhelio_rms_max_good, pyrhelio_rms_min::Real = pyrhelio_rms_min_good)
#	argmax(
#		df |> @map( sum((abs.(_.manifest.solar_hour_angle).<solar_hour_angle_threshold) .& map(i->_.manifest.order_snrs[i][order_idx_to_use_for_snr] > order_snr_threshold,1:size(_.manifest,1))) ) |> collect )
	argmax(
		df |> @map( sum((abs.(_.manifest.solar_hour_angle).<solar_hour_angle_threshold) .&
				( pyrhelio_ratio_min_good .<= _.manifest.pyrheliometer_flux_over_model.<= pyrhelio_ratio_max_good ) .&
				(pyrhelio_rms_min_good .<= _.manifest.pyrheliometer_rms_flux .<= pyrhelio_rms_max_good) )) 	|> collect)
end

# ╔═╡ 8cb9eca6-2ebe-4faa-b6c6-5ea267e64dc9
obs_idx_good_day = find_good_day_idx(datafiles_use) # 71 = 03/27

# ╔═╡ f2ed8d60-28e8-4fd4-8c08-88cf0e26ac1d
x = load(datafiles_use.filename[obs_idx_good_day]);

# ╔═╡ 8d6bda99-27e1-4993-9af5-816910d7bb21
# @time 
sigmarv3 = map(ordix->calc_ccf_order_σ_rv_expr2(x["v_grid"],x["order_ccfs"][:,ordix,100], x["order_ccf_vars"][:,ordix,100]),orders_plt)

# ╔═╡ 05d7024b-8a67-48c2-af78-527e9f1f8718
sqrt(1.0 ./ sum(1.0 ./ sigmarv3.^2))

# ╔═╡ 5a155b3a-c1b6-4edb-8c5d-14cd490f9ce8
function group_days(dates::AVD; max_days_in_bin::Integer = 7, min_obs_in_bin = 4 ) where { AVD<:AbstractVector{Date} }
	Δdays = Dates.value.(dates[2:end].-dates[1:end-1])
	num_obs_in_next_max_days = zeros(Int64,length(dates))
	idx_last_day_in_bin = zeros(Int64,length(dates))
	#start = 1
	#local stop 
	#groups = (Int64[])[]
	already_used = falses(length(dates))
	for i in 1:length(dates)
		mask_close = Dates.value.(dates[i:min(i+max_days_in_bin,length(dates))] .- dates[i]) .< max_days_in_bin
		mask = mask_close # .& .! already_used
		idx_last_day_in_bin[i] = findlast(mask)
		num_obs_in_next_max_days[i] = sum(mask)
	end
	possible_start_idx = findall(num_obs_in_next_max_days .>= min_obs_in_bin)
	possible_stop_idx =  possible_start_idx.+idx_last_day_in_bin[possible_start_idx].-1
	
	idx_keep = falses(length(possible_start_idx))
	
	for i in 1:length(possible_start_idx)
		if !already_used[possible_start_idx[i]] 
			idx_keep[i] = true
			for j in (possible_start_idx[i]+1):possible_stop_idx[i]
				already_used[j] = true
			end
		end
		
	end
	#return possible_start_idx[idx_keep], possible_stop_idx[idx_keep], num_obs_in_next_max_days[possible_start_idx[idx_keep]]
	return (idx_day_ranges = UnitRange.(possible_start_idx[idx_keep], possible_stop_idx[idx_keep]), 
		num_days = num_obs_in_next_max_days[possible_start_idx[idx_keep]] )
	#return possible_start_idx, possible_start_idx.+idx_last_day_in_bin[possible_start_idx].-1, num_obs_in_next_max_days[possible_start_idx]
	
	#return num_obs_in_next_max_days
	#return groups
end

# ╔═╡ 78f03589-85e8-402d-9474-4fe2e5b8caf4
day_idx_months, days_in_month = group_days(Dates.Date.(datafiles_use.datestr,"YYYYmmdd"),max_days_in_bin=30, min_obs_in_bin=9)

# ╔═╡ 49f1594c-c271-4981-9753-da88be755b23
day_idx_weeks, days_in_week = group_days(Dates.Date.(datafiles_use.datestr,"YYYYmmdd"),max_days_in_bin=7, min_obs_in_bin=4)

# ╔═╡ 019c4b03-caf9-48eb-b75d-882582a19d94
function calc_std_rv_within_period(dates::AVD, rvs::AVR; max_days_in_bin::Integer, min_obs_in_bin::Integer )  where { AVD<:AbstractVector{Date}, T1<:Real, AVR<:AbstractVector{T1} }
	@assert length(dates) == length(rvs)
	non_nans = findall(.!isnan.(rvs))
	day_idx_weeks, days_in_week = group_days(Dates.Date.(datafiles_use.datestr[non_nans],"YYYYmmdd"),max_days_in_bin=max_days_in_bin, min_obs_in_bin=min_obs_in_bin)
	num_weeks = length(day_idx_weeks)
	rmss_weekly = zeros(num_weeks)
	for i in 1:num_weeks
		rmss_weekly[i] = std(rvs[non_nans][day_idx_weeks[i]])
	end
	return rmss_weekly
end


# ╔═╡ 5140cedf-cdd4-43a3-8778-3f926a9f31c0
	function calc_weekly_std_rv(dates::AVD, rvs::AVR)  where { AVD<:AbstractVector{Date}, T1<:Real, AVR<:AbstractVector{T1} }
		calc_std_rv_within_period(dates,rvs, max_days_in_bin=7, min_obs_in_bin=4)
	end

# ╔═╡ c7690425-1a71-4fdf-9b41-e7f4e5f13d5b
	function calc_monthly_std_rv(dates::AVD, rvs::AVR)  where { AVD<:AbstractVector{Date}, T1<:Real, AVR<:AbstractVector{T1} }
		calc_std_rv_within_period(dates,rvs, max_days_in_bin=30, min_obs_in_bin=8)
	end

# ╔═╡ fe96f553-ec82-4997-a0ca-72f26f4aa5dd
md"# RESUME HERE"

# ╔═╡ 34bda3ea-4238-4670-9f6b-9207f63047cb
function calc_order_weights(df::DataFrame, obs_idx::Integer)
	order_weights = 1.0 ./ df.σ_rvs[obs_idx_good_day]
	order_weights[isinf.(order_weights)] .= 0
	order_weights[isnan.(order_weights)] .= 0
	order_weights = vec(sum(order_weights,dims=2))
	order_weights ./= sum(order_weights)
	return order_weights
end

# ╔═╡ 1174095b-a339-47b2-9f21-3f5b2c52b4eb
order_weights = calc_order_weights(all_rvs,obs_idx_good_day)

# ╔═╡ 72522f83-2b3c-4e40-a1de-c3508d50fdab
function summarize_day_orders(day_to_plt::Integer; df::DataFrame, orders, weights, solar_hour_angle_threshold::Real=solar_hour_angle_threshold_good, pyrhelio_ratio_min::Real = pyrhelio_ratio_min_good, pyrhelio_ratio_max::Real = pyrhelio_ratio_max_good, pyrhelio_rms_max::Real = pyrhelio_rms_max_good, pyrhelio_rms_min::Real = pyrhelio_rms_min_good )
	local shas = df.manifest[day_to_plt].solar_hour_angle ./24
	local times = df.manifest[day_to_plt].bjd
	local tsha0 = roots(Polynomials.fit(times.-2.459e6,shas,1))[1].+2.459e6

	local δrvs = df.manifest[day_to_plt][!,:Δv_diff_ext]
	sha_mask = abs.(df.manifest[day_to_plt].solar_hour_angle).<= solar_hour_angle_threshold
	pyrh_mask = (pyrhelio_ratio_min .<= df.manifest[day_to_plt].pyrheliometer_flux_over_model .<= pyrhelio_ratio_max ) .& ( pyrhelio_rms_min .<= df.manifest[day_to_plt].pyrheliometer_rms_flux .<= pyrhelio_rms_max )	
	mask = sha_mask .& pyrh_mask
	if sum(mask) ==0
		return Dict(:mean_order_rv=>fill(NaN,length(orders)),
		:rms_order_rv=>fill(NaN,length(orders)), :t_sha0=>tsha0, :Δt=>times.-tsha0, 
		:order_rv_slope=>fill(NaN,length(orders)), :order_rv_sha0=>fill(NaN,length(orders)),
		:rms_order_rv_binned=>fill(NaN,length(orders)) )
	end
	
	local daily_means = map(ord->NaNMath.mean(df.rvs[day_to_plt][ord,mask].+δrvs[mask]),orders)
	local daily_rmss = map(ord->NaNMath.std(df.rvs[day_to_plt][ord,mask].+δrvs[mask]),orders)	
	local daily_mean_σ = map(ord->NaNMath.mean(df.σ_rvs[day_to_plt][ord,mask]),orders)	
	daily_rmss_binned = fill(NaN,length(orders))
	for ord in orders
		(t_binned, rvs_binned, num_obs_in_bin) = bin_times_and_rvs_max_Δt(times=times[mask], 
			rvs=df.rvs[day_to_plt][ord,mask].+δrvs[mask], 
				Δt_threshold= 6.0/(24*60))
		mask_bin_full = num_obs_in_bin .== 4
		daily_rmss_binned[ord] = sum(mask_bin_full)>3 ? std(rvs_binned[mask_bin_full]) : NaN
	end
			
	local daily_fits = map(ord->Polynomials.fit(times[mask].-tsha0,df.rvs[day_to_plt][ord,mask].+δrvs[mask],1),orders)
	local daily_rv_sha0s = map(i->daily_fits[i][0],1:length(orders))
	local daily_slopes = map(i->daily_fits[i][1],1:length(orders))
	
	return Dict(:mean_order_rv=>daily_means, :mean_order_σrv=> daily_mean_σ,
		:rms_order_rv=>daily_rmss, :t_sha0=>tsha0, :Δt=>times.-tsha0, 
		:order_rv_slope=>daily_slopes, :order_rv_sha0=>daily_rv_sha0s,
		:rms_order_rv_binned=>daily_rmss_binned )
end

# ╔═╡ 4c62a111-112c-4b8d-979d-dea84e35232a
day_orders = DataFrame(summarize_day_orders.(1:size(datafiles_use,1),df=df_rvs,orders=1:num_order_idx, weights= order_weights,  pyrhelio_ratio_min=pyrhelio_ratio_min_good, pyrhelio_ratio_max=pyrhelio_ratio_max_good, pyrhelio_rms_max=pyrhelio_rms_max_good, pyrhelio_rms_min=pyrhelio_rms_min_good) )

# ╔═╡ 2d56d9c0-b23a-4399-a86e-77edbf25d8b5
median_order_rv = map(ord->NaNMath.median(map(i->day_orders.mean_order_rv[i][ord],(1:size(day_orders,1)))),1:num_order_idx)

# ╔═╡ 3d088149-daa9-4664-9a06-7685d4cbb00e
day_orders.mean_order_rv[obs_idx_good_day]

# ╔═╡ 701ebc68-2809-4544-be85-c4fadae6c391
begin
	daily_rms_across_orders = map(obsid->NaNMath.std((day_orders.order_rv_sha0[obsid].-median_order_rv)[vcat(1,4:6,8:10,12)]),1:size(day_orders,1))
end

# ╔═╡ ea26ca76-102d-4011-9f8b-0917c021cc9f
all_rvs[5,:rvs]

# ╔═╡ cf17cae1-69ea-4050-b302-cfc5d30d6235
md"RESUME 2" 

# ╔═╡ 6fb4b458-6192-4000-8efd-f6fb02488861
function summarize_day_total(day_to_plt::Integer; df::DataFrame, orders, weights, 
		solar_hour_angle_threshold::Real=solar_hour_angle_threshold_good, pyrhelio_ratio_min::Real = pyrhelio_ratio_min_good, pyrhelio_ratio_max::Real = pyrhelio_ratio_max_good, pyrhelio_rms_max::Real = pyrhelio_rms_max_good, pyrhelio_rms_min::Real = pyrhelio_rms_min_good)
	#orders_to_use = findall(vec(sum(isnan.(all_rvs.rvs[day_to_plt][orders,:]),dims=2).==0))
	
	#obs_idx_good = vec(sum(isnan.(all_rvs[day_to_plt,:rvs]),dims=1)).<1
	
	rv_tot = vec(sum(df.rvs[day_to_plt][orders,:].*weights[orders],dims=1) ./ 					   sum(weights[orders],dims=1))
	σ_rv_tot = vec(sum(df.σ_rvs[day_to_plt][orders,:].*weights[orders],dims=1) ./ 					   sum(weights[orders],dims=1))  #.* sqrt(NEID.default_ccf_mask_v_width(NEID2D())/(v_grid[2]-v_grid[1]))
	rv_tot .+= df.manifest[day_to_plt][:,:Δv_diff_ext]
	times = df.manifest[day_to_plt].bjd # [obs_idx_good]
	shas = df.manifest[day_to_plt].solar_hour_angle #=[obs_idx_good]=# ./24
	tsha0 = roots(Polynomials.fit(times.-2.459e6,shas,1))[1].+2.459e6

	sha_mask = abs.(df.manifest[day_to_plt].solar_hour_angle #=[obs_idx_good]=# ).<= solar_hour_angle_threshold
	#snr_mask = map(i->df.manifest[day_to_plt].order_snrs[i][order_idx_to_use_for_snr],1:length(df.manifest[day_to_plt].order_snrs)) .>= order_snr_threshold
	pyrh_mask = (pyrhelio_ratio_min .<= df.manifest[day_to_plt].pyrheliometer_flux_over_model .<= pyrhelio_ratio_max ) .& ( pyrhelio_rms_min .<= df.manifest[day_to_plt].pyrheliometer_rms_flux .<= pyrhelio_rms_max )	
	#mask = sha_mask .& snr_mask
	mask = sha_mask .& pyrh_mask
	if sum(mask) >= 1
		rms_rv = std(rv_tot[mask])
		(t_binned, rvs_binned, num_obs_in_bin) = 
			bin_times_and_rvs_max_Δt(times=times[mask],rvs=rv_tot[mask],Δt_threshold= 6.01/(24*60))
		mask_bin_full = num_obs_in_bin .== 4
		rms_rv_binned = sum(mask_bin_full)>3 ? std(rvs_binned[mask_bin_full]) : NaN
		daily_fit = Polynomials.fit(times[mask].-tsha0,rv_tot[mask],1)
		rms_rv_binned_fit =  sum(mask_bin_full)>3 ? std(rvs_binned[mask_bin_full] .-daily_fit.(t_binned[mask_bin_full].-tsha0)) : NaN
		fit_constant = daily_fit[0]
		fit_slope = daily_fit[1]
		num_obs_in_bin = zeros(Int64,0)
	else
		rms_rv = NaN
		t_binned = zeros(0)
		rvs_binned = zeros(0)
		num_obs_in_bin = zeros(Int64,0)
		rms_rv_binned = NaN
		rms_rv_binned_fit = NaN
		fit_constant = NaN
		fit_slope = NaN
	end

	Dict( :rv=>rv_tot, :σ_rv=>σ_rv_tot, :rms_rv=>rms_rv,
		  :t_sha0=>tsha0, :Δt => times.-tsha0, 
		  :rv_slope => fit_slope, :rv_sha0 =>fit_constant,
		  :Δt_binned => t_binned.-tsha0, :rv_binned=>rvs_binned, :num_obs_in_bin=>num_obs_in_bin,
		  :rms_rv_binned=>rms_rv_binned, :rms_rv_binned_fit =>rms_rv_binned_fit,
		  :mask=>mask, :wavelength_cal_flag => df.lfc_flag[day_to_plt] || df.drift_flag[day_to_plt]
		)
end

# ╔═╡ 343f3226-a6e1-4bf5-a816-f92cb74cf434
begin
	day_comb = DataFrame(summarize_day_total.(1:size(datafiles_use,1),df=df_rvs,orders=order_idx_for_rvs, weights= order_weights, pyrhelio_ratio_min=pyrhelio_ratio_min_good, pyrhelio_ratio_max=pyrhelio_ratio_max_good, pyrhelio_rms_max=pyrhelio_rms_max_good, pyrhelio_rms_min=pyrhelio_rms_min_good));
end;

# ╔═╡ 706da46c-044c-49e5-9a77-103afe913797
begin
	local plt = plot(legend=:none)
	histogram!(plt,map(i->size(day_comb.rv_binned[i],1),obs_idx_all_good_days),nbins=40)
	xlabel!(plt,"Number of binned observations")
	ylabel!(plt,"Number of days")
	title!(plt,"Distribution of binned observations/day")
end

# ╔═╡ 2fd9c2e4-fd39-44b7-8d88-bddaf46ee26b
day_comb.σ_rv

# ╔═╡ 592b8b00-91a0-4ab4-b26b-9fc358978e9a
md"""
Day to plot
$(@bind obs_idx_plt NumberField(1:size(day_comb,1); default=obs_idx_good_day))  
"""

# ╔═╡ 72ad7501-5a9c-4e06-91e9-a03baff07a83
size.(day_comb.rv_binned[obs_idx_all_good_days],1)

# ╔═╡ a782e454-470f-4373-ad93-2b73c0a69b8b
neid_7day = calc_weekly_std_rv(Dates.Date.(datafiles_use.datestr,"YYYYmmdd"),day_comb.rv_sha0)

# ╔═╡ 6a5425e0-7fbf-4e62-b839-bb34f00bd026
begin
	local plt = plot()
	local bins = range(0,1.6,length=20)
	histogram!(plt,harpsn_7day.rms_rv,bins=bins,alpha=0.5,label="HARPS-N")
	histogram!(plt,neid_7day, alpha=0.5, bins=bins,label="NEID")
	xlabel!(plt,"RMS RV (m/s)")
	ylabel!(plt,"Count")
	title!(plt,"RMS RV over 7 day periods")
	if save_figs  savefig("rms_rvs_7day_histo.png") end
	plt
end

# ╔═╡ 7aa5561f-682a-4c15-b6d3-4c3ef8572da1
neid_30day = calc_monthly_std_rv(Dates.Date.(datafiles_use.datestr,"YYYYmmdd"),day_comb.rv_sha0)

# ╔═╡ a46da6b4-c8c9-4c1f-9bb6-33899bc63336
begin
	local plt = plot()
	local bins = range(0,2.0,length=20)
	histogram!(plt,harpsn_30day.rms_rv,bins=bins,alpha=0.5,label="HARPS-N")
	histogram!(plt,neid_30day, alpha=0.5, bins=bins,label="NEID")
	xlabel!(plt,"RMS RV (m/s)")
	ylabel!(plt,"Count")
	title!(plt,"RMS RV over 30 day periods")
	if save_figs  savefig("rms_rvs_30day_histo.png") end
	plt
end

# ╔═╡ c863a4fc-270b-4473-beaa-d34ba2d719de
begin
	mask_decent_days = (length.(day_comb.rv_binned).>=min_obs_binned_in_day) .& .! day_comb.wavelength_cal_flag # .& (daily_rms_across_orders .< 6) #.&
	all_days_to_plt = plt_all_days ? (1:size(day_orders,1)) : (1:size(day_orders,1))[mask_decent_days]
	num_all_days_to_plt = length(all_days_to_plt)
end

# ╔═╡ fb17d03b-880c-4be6-bd50-f31bdd3e2fd4
md"""
Date range (min, max, step) 
$(@bind days_to_plt_idx_lo Slider(1:num_all_days_to_plt; default=1))  
$(@bind days_to_plt_idx_hi Slider(1:num_all_days_to_plt; default=num_all_days_to_plt))

"""

# ╔═╡ bbf59177-36a8-480c-af39-1998b37032d2
begin
		default_order_rv_lo = plt_all_days ? -500 : -6
		default_order_rv_hi = plt_all_days ? 500 : 6
		default_order_rv_step = plt_all_days ? 1 : 0.1
		days_to_plt = all_days_to_plt[days_to_plt_idx_lo:days_to_plt_idx_hi]
end;

# ╔═╡ 92175df7-b060-4f6f-8acb-265bdbe6713a
begin
	max_order_daily_rms_hi = NaNMath.maximum(daily_rms_across_orders[days_to_plt])
	default_order_daily_rms_hi = 10 #NaNMath.maximum(daily_rms_across_orders[mask_decent_days])
end;

# ╔═╡ 8af9d40f-8d50-4376-82e4-01981a54893f
md"""
Order RV RMS limits (lo, hi) 
$(@bind order_daily_rms_lo Slider(0.0:0.2:5; default=0.1))  
$(@bind order_daily_rms_hi Slider(0.0:0.2:30; default=default_order_daily_rms_hi))

"""

# ╔═╡ f7db71f7-9bac-4aad-a636-5c60ae4c259b
begin
	local plt = plot()
	local pal = palette(:berlin, length(orders_plt))
	for (i,ordidx) in enumerate(orders_plt)
		scatter!(plt,Date.(Dates.julian2datetime.(day_orders.t_sha0[days_to_plt])),	map(obsidx->day_orders[obsidx,:rms_order_rv][ordidx],1:size(day_orders,1)),ms=1.0,label=:none,c=pal[i])
		scatter!(plt,Date.(Dates.julian2datetime.(day_orders.t_sha0[days_to_plt])),	map(obsidx->day_orders[obsidx,:rms_order_rv][ordidx],days_to_plt),label=:none,c=pal[i])
	end
	xlabel!(plt,"Time (d)")
	ylabel!(plt,"RMS of Orders RVs within each day")
	ylims!(order_daily_rms_lo,order_daily_rms_hi)
	title!(plt,"Daily RMS of Order RVs")	
end

# ╔═╡ f57e8544-204a-47a0-93c9-24a2982a0aed
begin
	local plt = plot(legend=:none)
	ylims!(plt,order_daily_rms_lo,order_daily_rms_hi)
	for i in 1:length(day_orders.rms_order_rv)
		scatter!(plt,174 .-(first(orders_to_use_default(NEID2D())).+orders_plt.-1),day_orders.rms_order_rv[i][orders_plt], label=string(i))
	end
	#scatter!(orders_plt,day_orders.mean_order_σrv[1][orders_plt], label="Mean σ")
	#ylims!(0,10)
	xlabel!(plt,"Physical Order")
	ylabel!(plt,"RMS (m/s)")
	title!(plt,"RMS of Daily RVs versus Order")
	plt
end

# ╔═╡ f8e4b341-1613-4993-a193-437f2adf4a21
md"""
Order RV limits (lo, hi) $(@bind order_rv_lo Slider(-200.0:0.0; default=default_order_rv_lo))  
$(@bind order_rv_hi Slider(0.0:200.0; default=default_order_rv_hi))
"""

# ╔═╡ 0d9dd1ea-9ef8-460c-8d95-20a46eb1ebb6
begin
	local plt1 = plot(legend=:none)
	local t_day =  map(i->day_orders[i,:t_sha0],1:size(day_orders,1))
	local orders_to_plt =use_order_range ?  (order_idx_lo:order_idx_step:order_idx_hi) :  order_idx_for_rvs
	local pal = palette(:berlin, length(orders_to_plt))
	for (i,ord) in enumerate(orders_to_plt)
		ord_rv_vs_day = map(i->day_orders[i,:order_rv_sha0][ord]-median_order_rv[ord],days_to_plt)
		if all(ord_rv_vs_day .== 0) continue end
		#scatter!(plt1,Date.(Dates.julian2datetime.(t_day)),map(i->day_orders[i,:order_rv_sha0][ord]-median_order_rv[ord],1:size(day_orders,1)),c=pal[i],ms=1.0)
		scatter!(plt1,Date.(Dates.julian2datetime.(t_day[days_to_plt])),ord_rv_vs_day,c=pal[i])
	end
	title!(plt1,"Order RV at midday")
	xlabel!(plt1,"BJD")
	ylabel!(plt1,"ΔRV (m/s)")
	ylims!(plt1,order_rv_lo,order_rv_hi)
	plt1
end

# ╔═╡ cfce959c-c171-4923-b5f9-29666d9c20e4
begin
	local plt2 = plot(legend=:none)
	local t_day =  map(i->day_orders[i,:t_sha0],1:size(day_orders,1))
	local orders_to_plt = use_order_range ? (order_idx_lo:order_idx_step:order_idx_hi) : order_idx_for_rvs
	local pal = palette(:berlin, length(orders_to_plt))
	for (j,ord) in enumerate(orders_to_plt)
		ord_rv_vs_day = map(i->day_orders[i,:order_rv_slope][ord],days_to_plt)
		if all(ord_rv_vs_day .== 0) continue end
		scatter!(plt2,Date.(Dates.julian2datetime.(t_day[days_to_plt])),ord_rv_vs_day,c=pal[j])
	end
	title!(plt2,"Order RV Daily Slope")
	xlabel!(plt2,"BJD")
	ylabel!(plt2,"RV (m/s/day)")
	ylims!(plt2,order_slope_lo,order_slope_hi)
	plt2
end

# ╔═╡ 0b07ac3c-d61e-43ff-acff-f907019835d6
begin
	local plt = plot(legendtitle="Order Index")
	local pal = palette(:berlin, length(order_idx_for_rvs))
	for (i,ord) in enumerate(order_idx_for_rvs)	
		#histogram!(plt,day_orders[ord,:rms_order_rv_binned],bin=range(0.1,stop=order_rv_rms_hi,length=40), color=pal[i],label=string(ord))
		histogram!(plt,map(obsid->day_orders[obsid,:rms_order_rv_binned][ord],days_to_plt),
				bin=range(0.1,stop=order_rv_rms_hi,length=40), color=pal[i],alpha=0.4,label=string(ord))
	end
	if 1<= order_to_highlight <= size(day_orders,1)
		histogram!(plt,map(obsid->day_orders[obsid,:rms_order_rv_binned][order_to_highlight],days_to_plt),
				bin=range(0.1,stop=order_rv_rms_hi,length=40), color=:yellow,alpha=0.5,label=string(order_to_highlight))
	end
	xlabel!(plt,"RMS (m/s)")
	ylabel!(plt,"Count")
	title!(plt,"RMS of order RVs (6min bin) within each good day")
	xlims!(plt,order_rv_rms_lo,order_rv_rms_hi)
	if save_figs savefig(plt,"rms_order_rv_binned.png") end
	plt
end

# ╔═╡ 9dc77012-f845-4a5b-a059-fd13c1308b3d
begin
	local mask = .! day_comb.wavelength_cal_flag
	
	local plt = plot(legend=:none) 
#	histogram!(plt,daily_rms_across_orders[mask], 
#		bins=range(order_daily_rms_lo,stop=order_daily_rms_hi,length=20),alpha=0.3, label="All usable days")
	histogram!(plt,daily_rms_across_orders[days_to_plt], 
		bins=range(order_daily_rms_lo,stop=order_daily_rms_hi,length=20),alpha=0.3, label="Good days")
	xlabel!(plt,"RMS (m/s)")
	ylabel!(plt,"Count")
	title!(plt,"RMS in daily RV over orders used for RVs")
	xlims!(plt,order_daily_rms_lo,order_daily_rms_hi)
	plt
end

# ╔═╡ e858ca22-5ddd-45aa-8815-c3acc3f9741e
date_mask = datetime2julian(Date(2021,01,01)) .<= day_comb.t_sha0[days_to_plt] .<= datetime2julian(Date(2021,05,15));

# ╔═╡ c8621274-00c9-4f67-b285-c624d1346527
begin
	default_daily_rv_rms_lo = NaNMath.minimum(day_comb.rms_rv_binned_fit)
	default_daily_rv_rms_hi = NaNMath.maximum(day_comb.rms_rv)
end;

# ╔═╡ 9b81f986-ea3c-49f5-be04-b12b392aacd6
md"""
Daily RMS RV limits (lo, hi)
$(@bind daily_rv_rms_lo Slider(0.0:0.05:1.0; default=default_daily_rv_rms_lo))
$(@bind daily_rv_rms_hi Slider(0.0:0.2:5.0; default=1.2))
"""

# ╔═╡ 34397405-e939-464c-a384-5d6c8fb9ca83
begin
	local plt = histogram(day_comb.rms_rv[days_to_plt], bins=daily_rv_rms_lo:(round((daily_rv_rms_hi-daily_rv_rms_lo)/40,digits=3)):daily_rv_rms_hi, alpha=1.0,label="Unbinned")
	histogram!(plt,day_comb.rms_rv_binned[days_to_plt], bins=daily_rv_rms_lo:(round((daily_rv_rms_hi-daily_rv_rms_lo)/40,digits=3)):daily_rv_rms_hi,alpha=1.0, color=:red, label="Binned (6 min)")
	#histogram!(day_comb.rms_rv_binned_fit[days_to_plt], bins=daily_rv_rms_lo:(round((daily_rv_rms_hi-daily_rv_rms_lo)/40,digits=3)):daily_rv_rms_hi, alpha=0.4, label="Binned - linear fit")
	
	xlims!(plt,daily_rv_rms_lo,daily_rv_rms_hi)
	xlabel!(plt,"RMS (m/s)")
	ylabel!(plt,"Count")
	title!(plt,"RMS of weighted RVs within each good day")
	if save_figs savefig(plt,"rms_rv_binned.png") end
	plt
end

# ╔═╡ 7b2e6266-6dfe-47b5-8f28-3b9d76f57a42
begin
	plt = plot(legend=:none)
	scatter!(plt,Date.(Dates.julian2datetime.(day_comb.t_sha0)),day_comb.rms_rv,ms=1.5,color=:blue)
	scatter!(plt,Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt])),day_comb.rms_rv[days_to_plt], color=:blue)
	scatter!(plt,Date.(Dates.julian2datetime.(day_comb.t_sha0)),day_comb.rms_rv_binned, ms=2,color=:red)
	scatter!(plt,Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt])),day_comb.rms_rv_binned[days_to_plt], color=:red)
	title!(plt,"Daily RV RMS")
	xlabel!(plt,"Day")
	ylabel!(plt,"RMS of weighted RVs within each good day (m/s)")
	xlims!(extrema(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt]))))
	ylims!(plt,daily_rv_rms_lo,daily_rv_rms_hi)
	if save_figs savefig(plt,"rms_rv_vs_time.png") end
	plt
end

# ╔═╡ 7d63e3eb-6132-439a-bc68-05e73cd5b068
Δrv = day_comb.rv_sha0.-NaNMath.median(day_comb.rv_sha0[days_to_plt][date_mask])	

# ╔═╡ d4e7987c-94bf-4194-8e16-88c2581832d1
begin
	default_daily_rv_lo = NaNMath.minimum(Δrv[days_to_plt])
	default_daily_rv_hi = NaNMath.maximum(Δrv[days_to_plt])
end;

# ╔═╡ 72e90532-02da-4d31-88c3-7766393d9f92
md"""
Daily RV limits (lo, hi)
$(@bind daily_rv_lo Slider(-20.0:0.2:0.0; default=default_daily_rv_lo))
$(@bind daily_rv_hi Slider(0.0:0.2:20.0; default=default_daily_rv_hi))

"""

# ╔═╡ fca2fa7d-99fe-44bd-a34a-1b1552a3ae48
begin
	local plt = plot(legend=:none)
	local fit_mask = date_mask .&  (-10.0 .<Δrv[days_to_plt].< 10.0 )
	long_term_daily_rv_fit = Polynomials.fit(day_comb.t_sha0[days_to_plt][fit_mask],Δrv[days_to_plt][fit_mask],1)
	scatter!(plt,Date.(Dates.julian2datetime.(day_comb.t_sha0)),Δrv,ms=2,color=:black)
	scatter!(plt,Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt][fit_mask])), Δrv[days_to_plt])
	local pred = long_term_daily_rv_fit.(day_comb.t_sha0)
	#plot!(plt,Date.(Dates.julian2datetime.(day_comb.t_sha0)),pred)
	rms_resid = std(Δrv[days_to_plt][fit_mask].-pred[days_to_plt][fit_mask])
	rms_resid_0 = std(Δrv[days_to_plt][fit_mask])
	title!(plt,"RV fit at solar noon.  " * 
		"RMS: " * string(round(rms_resid_0,digits=3)) * "m/s" # *
		#"RMS (fit): " * string(round(rms_resid,digits=3)) * "m/s"  * 
		#"\nBest-fit slope: " * string(round(long_term_daily_rv_fit[1]*100,digits=2)) * " cm/s/day"
	)
	xlabel!(plt,"Day")
	ylabel!(plt,"RV at solar noon (m/s)")
	xlims!((minimum(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt].-1))),maximum(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt].+1)))))
	ylims!(plt,daily_rv_lo,daily_rv_hi)
	#title!(plt,"RV Fit at Solar Noon")
	if save_figs savefig(plt,"rv_sha0_vs_time.png") end
	plt
end

# ╔═╡ 1c42a7ce-3e97-48b3-8ddd-89cab6177bfe
"Daily RVs:  Linear slope: " * string(round(long_term_daily_rv_fit[1],digits=4)) * " m/s/d    " * string("RMS about fit: " * string(round(rms_resid,digits=3)) * " m/s") * " Median RMS Binned RVs " * string(round(NaNMath.median(day_comb.rms_rv_binned[days_to_plt]),digits=3)) * " m/s"

# ╔═╡ ea5aa2da-f64f-44a7-8c41-746c9e0145e7
"Daily RVs:  Linear slope: " * string(round(long_term_daily_rv_fit[1],digits=4)) * " m/s/d    " * string("RMS about fit: " * string(round(rms_resid,digits=3)) * " m/s") * " Median RMS Binned RVs " * string(round(NaNMath.median(day_comb.rms_rv_binned[days_to_plt]),digits=3)) * " m/s"

# ╔═╡ 0e9b78b0-208f-4b29-96c3-4b4129416bdd
md"""
Daily RV slope limits (lo, hi)
$(@bind daily_slope_lo Slider(-20.0:0.0; default=minimum(day_comb.rv_slope[days_to_plt])))
$(@bind daily_slope_hi Slider(0.0:20.0; default=maximum(day_comb.rv_slope[days_to_plt])))

"""

# ╔═╡ 4c71c3d6-1f9e-4812-9f36-2e65a43c8207
begin
	plot(legend=:none)
	scatter!(Date.(Dates.julian2datetime.(day_comb.t_sha0)),day_comb.rv_slope,ms=2,color=:black)
	scatter!(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt])),day_comb.rv_slope[days_to_plt])
	title!("Daily RV slope")
	xlabel!("Day")
	ylabel!("RV Slope (m/s/day)")
	xlims!(extrema(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt]))))
	ylims!(daily_slope_lo,daily_slope_hi)
end

# ╔═╡ b2b81585-2144-4b37-b0a7-286c862bb812
begin
	plot(legend=:none)
	#histogram(day_comb.rv_slope[days_to_plt],bins=range(daily_slope_lo,stop=daily_slope_hi,length=30),alpha=0.5,label=:none)
	histogram(day_comb.rv_slope[days_to_plt],nbins=30,alpha=0.5,label=:none)
	#local idx_good = length.(day_comb.rv_binned[days_to_plt]).>=25
	#histogram!(day_comb.rv_slope[days_to_plt][idx_good],bins=range(daily_slope_lo,stop=daily_slope_hi,length=30),alpha=0.5,label="Many points")
	title!("Distribution of Daily RV slopes")
	xlabel!("RV Slope (m/s/day)")
	xlims!(daily_slope_lo,daily_slope_hi)
end

# ╔═╡ e9f30692-375c-49fc-bf0d-f045951ae11f
begin
	#=
	# DRP 0.7
	day_farred = DataFrame(summarize_day_total.(1:size(datafiles_use,1),df=df_rvs,orders=vcat(37:41,46:50), weights= order_weights,  pyrhelio_ratio_min=pyrhelio_ratio_min_good, pyrhelio_ratio_max=pyrhelio_ratio_max_good, pyrhelio_rms_max=pyrhelio_rms_max_good, pyrhelio_rms_min=pyrhelio_rms_min_good))
	day_red = DataFrame(summarize_day_total.(1:size(datafiles_use,1),df=df_rvs,orders=vcat(13:14,16:21,23:25), weights= order_weights,  pyrhelio_ratio_min=pyrhelio_ratio_min_good, pyrhelio_ratio_max=pyrhelio_ratio_max_good, pyrhelio_rms_max=pyrhelio_rms_max_good, pyrhelio_rms_min=pyrhelio_rms_min_good))	
	day_blue = DataFrame(summarize_day_total.(1:size(datafiles_use,1),df=df_rvs,orders=vcat(1:6,8,10:12), weights= order_weights,  pyrhelio_ratio_min=pyrhelio_ratio_min_good, pyrhelio_ratio_max=pyrhelio_ratio_max_good, pyrhelio_rms_max=pyrhelio_rms_max_good, pyrhelio_rms_min=pyrhelio_rms_min_good))
	=#
	# DRP 1.0
	day_farred = DataFrame(summarize_day_total.(1:size(datafiles_use,1),df=df_rvs,orders=vcat(37:43,45:50), weights= order_weights,  pyrhelio_ratio_min=pyrhelio_ratio_min_good, pyrhelio_ratio_max=pyrhelio_ratio_max_good, pyrhelio_rms_max=pyrhelio_rms_max_good, pyrhelio_rms_min=pyrhelio_rms_min_good))
	day_red = DataFrame(summarize_day_total.(1:size(datafiles_use,1),df=df_rvs,orders=vcat(18:21,23:29), weights= order_weights,  pyrhelio_ratio_min=pyrhelio_ratio_min_good, pyrhelio_ratio_max=pyrhelio_ratio_max_good, pyrhelio_rms_max=pyrhelio_rms_max_good, pyrhelio_rms_min=pyrhelio_rms_min_good))	
	day_blue = DataFrame(summarize_day_total.(1:size(datafiles_use,1),df=df_rvs,orders=vcat(1:15), weights= order_weights,  pyrhelio_ratio_min=pyrhelio_ratio_min_good, pyrhelio_ratio_max=pyrhelio_ratio_max_good, pyrhelio_rms_max=pyrhelio_rms_max_good, pyrhelio_rms_min=pyrhelio_rms_min_good))
	
end;

# ╔═╡ e18232e8-31e6-4fbf-8c96-5c36c9e8602b
begin
	Δrv_blue = day_blue.rv_sha0.-NaNMath.mean(day_blue.rv_sha0[days_to_plt])	
	Δrv_red =  day_red.rv_sha0.-NaNMath.mean(day_red.rv_sha0[days_to_plt])	
	Δrv_farred =  day_farred.rv_sha0.-NaNMath.mean(day_farred.rv_sha0[days_to_plt])	
end;

# ╔═╡ ab2e501d-1c2d-48e1-be40-942cb4e66bc8
md"""
Daily ΔRV limits (lo, hi)
$(@bind daily_Δrv_lo Slider(-20.0:0.0; default=minimum(vcat(Δrv_red[days_to_plt].-Δrv_blue[days_to_plt],Δrv_farred[days_to_plt].-Δrv_blue[days_to_plt]))-1))
$(@bind daily_Δrv_hi Slider(0.0:20.0; default=maximum(vcat(Δrv_red[days_to_plt].-Δrv_blue[days_to_plt],Δrv_farred[days_to_plt].-Δrv_blue[days_to_plt]))+1))

"""

# ╔═╡ 9d89d84d-f107-4a8d-be18-55af8edb546e
begin
	local plt = plot(legend=:bottomright)
	xlims!(extrema(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt]))))
	ylims!(plt,daily_Δrv_lo,daily_Δrv_hi)
	
	#scatter!(plt,Date.(Dates.julian2datetime.(day_blue.t_sha0[days_to_plt])), Δrv_blue[days_to_plt],color=:blue)
	scatter!(plt,Date.(Dates.julian2datetime.(day_red.t_sha0[days_to_plt])), Δrv_red[days_to_plt].-Δrv_blue[days_to_plt],color=:green, label="Red-Blue")
	scatter!(plt,Date.(Dates.julian2datetime.(day_farred.t_sha0[days_to_plt])), Δrv_farred[days_to_plt].-Δrv_blue[days_to_plt],color=:red, label="FarRed-Blue")
	local pfit = Polynomials.fit(day_comb.t_sha0[days_to_plt],Δrv[days_to_plt],1)
	#rms_resid = std(Δrv[days_to_plt].-pfit.(day_comb.t_sha0[days_to_plt]))
	title!(plt,"ΔRV fit at midday") # \nLinear slope: " * string(round(pfit[1],digits=4)) * " m/s/d.  RMS about fit: " * string(round(rms_resid,digits=3)) * "m/s")
	xlabel!(plt,"Day")
	ylabel!(plt,"ΔRV at solar noon (m/s)")
	if save_figs savefig("daily_chromatic_rvs.png") end
	plt
end

# ╔═╡ 04e04706-cf7c-428d-8e3e-70643cf7ed60
md"""
Daily RV slope limits (lo, hi)
$(@bind daily_slope_subsets_lo Slider(-20.0:0.0; default=minimum(vcat(day_farred.rv_slope[days_to_plt],day_blue.rv_slope[days_to_plt],day_red.rv_slope[days_to_plt]))))
$(@bind daily_slope_subsets_hi Slider(0.0:20.0; default=maximum(vcat(day_farred.rv_slope[days_to_plt],day_blue.rv_slope[days_to_plt],day_red.rv_slope[days_to_plt]))))

"""

# ╔═╡ 833476e3-f1f1-4042-87e0-4c21ef37f7a7
begin
	plot(legend=:none)
	scatter!(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt])),day_comb.rv_slope[days_to_plt],ms=2,color=:black)
	scatter!(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt])),day_blue.rv_slope[days_to_plt],color=:blue)
	scatter!(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt])),day_red.rv_slope[days_to_plt],color=:green)
	scatter!(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt])),day_farred.rv_slope[days_to_plt],color=:red)
	title!("Daily RV slope for Order subsets")
	xlabel!("Day")
	ylabel!("RV Slope (m/s/day)")
	xlims!(extrema(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt]))))
	ylims!(daily_slope_subsets_lo,daily_slope_subsets_hi)
end

# ╔═╡ 31f212b1-dac3-4c77-a340-dd2837fc8807
function plot_rvs_day(obs_idx::Integer; solar_hour_angle_threshold::Real=solar_hour_angle_threshold_good, pyrhelio_ratio_min::Real = pyrhelio_ratio_min_good, pyrhelio_ratio_max::Real = pyrhelio_ratio_max_good, pyrhelio_rms_max::Real = pyrhelio_rms_max_good, pyrhelio_rms_min::Real = pyrhelio_rms_min_good)
	local plt = plot()
	sha_mask = abs.(datafiles_use.manifest[obs_idx].solar_hour_angle).<= solar_hour_angle_threshold
	 #snr_mask = map(i->datafiles_use.manifest[obs_idx].order_snrs[i][order_idx_to_use_for_snr],1:length(datafiles_use.manifest[obs_idx].order_snrs)) .>= order_snr_threshold
	pyrh_mask = (pyrhelio_ratio_min .<= datafiles_use.manifest[obs_idx].pyrheliometer_flux_over_model .<= pyrhelio_ratio_max ) .& ( pyrhelio_rms_min .<= datafiles_use.manifest[obs_idx].pyrheliometer_rms_flux .<= pyrhelio_rms_max )	
	#mask = sha_mask .& snr_mask
	mask = sha_mask .& pyrh_mask
	ytmp = day_comb.Δt[obs_idx][mask].*day_comb.rv_slope[obs_idx].+day_comb.rv_sha0[obs_idx]
	if sum(.!isnan.(ytmp)) >=3
		plot!(plt,day_comb.Δt[obs_idx][mask].*24,day_comb.Δt[obs_idx][mask].*day_comb.rv_slope[obs_idx].+day_comb.rv_sha0[obs_idx], label=:none)
	end
	if sum(.!isnan.(day_comb.rv[obs_idx][mask]))>0
		scatter!(plt,day_comb.Δt[obs_idx][mask].*24,	day_comb.rv[obs_idx][mask], yerr=day_comb.σ_rv[obs_idx], ms=1.5,color=:red,label="Unbinned")
		
	end	
	if sum(mask) >= 3
		(t_binned, rvs_binned, num_in_bin) = bin_times_and_rvs_max_Δt(times=day_comb.Δt[obs_idx][mask],rvs=day_comb.rv[obs_idx][mask],Δt_threshold= 5/(24*60))
		plt_mask = num_in_bin .== 4
		if sum(.!isnan.(rvs_binned[plt_mask])) >=1
			scatter!(plt,t_binned[plt_mask].*24,	rvs_binned[plt_mask], ms=5.0, color=:green, label="Binned 6min")
		end
	end
	warn_str = ""
	if datafiles_use.lfc_flag[obs_idx] 
		warn_str = warn_str * " LFC Flag"
	end
	if  datafiles_use.drift_flag[obs_idx]
		warn_str = warn_str * " Drift Flag"
	end
	xlabel!(plt,"Δt from solar noon (hr)")
	ylabel!(plt,"RV (m/s)" * warn_str)
	xlims!(plt,-solar_hour_angle_threshold,solar_hour_angle_threshold)
	#title!(plt,datafiles_use.datestr[obs_idx] * ": RMS= " * string(round(day_comb.rms_rv_binned[obs_idx],digits=3)) * "m/s Slope= " * string(round(day_comb.rv_slope[obs_idx],digits=3)) * "m/s/d" )
	return plt
end

# ╔═╡ 8bb918eb-cdb6-4a15-a344-88ec22df92b2
begin
	local plt = plot_rvs_day(obs_idx_plt, solar_hour_angle_threshold=solar_hour_angle_threshold_good, pyrhelio_ratio_min = pyrhelio_ratio_min_good, pyrhelio_ratio_max = pyrhelio_ratio_max_good, pyrhelio_rms_max = pyrhelio_rms_max_good, pyrhelio_rms_min = pyrhelio_rms_min_good)
	title!(plt,"Solar RVs " * datafiles_use.datestr[obs_idx_plt])
	xlims!(plt,-solar_hour_angle_threshold_good-0.5,solar_hour_angle_threshold_good+0.5)
	plt
end

# ╔═╡ 4441d217-5f2e-4052-8fda-e11672ae1813
md"""
## Package Setup
"""

# ╔═╡ e0756a3b-998a-4f8d-9ab3-df50422ce14b
md"""
## UI Tweaks
"""

# ╔═╡ 44ae02ac-6f9a-4af5-8901-3428da03af22
PlutoUI.TableOfContents()

# ╔═╡ b05ae27f-13df-4ae6-8189-1c21d86ee2df
html"""<style>
main {
    max-width: 900px;
    align-self: flex-start;
    margin-left: 20px;
}
"""

# ╔═╡ e3b8c1b9-cdb6-4a06-bc39-180f4bc42318


# ╔═╡ Cell order:
# ╟─87597e02-19c7-418a-8eb3-ad25695bb5bf
# ╟─ba3a32ee-9664-4279-b27b-0fed992e4d1c
# ╟─f64f9468-25ee-4799-b776-b1cfbe5b7f57
# ╠═c259acbf-da3d-4b73-af2b-8fc22f2548d7
# ╠═15741d2d-5fc2-4dc0-baf1-3ba3a72b5399
# ╟─d6219a57-7c1d-4b23-9b80-ae96f005b762
# ╠═dfea1ed2-4431-4d5b-99c7-560dc4472ab1
# ╠═011c5463-2edf-4903-b512-0df12d8b3779
# ╟─7f12ff80-e32c-400b-8cdb-1a800fdbc791
# ╟─1c42a7ce-3e97-48b3-8ddd-89cab6177bfe
# ╟─706da46c-044c-49e5-9a77-103afe913797
# ╟─1ee5e5d4-9b2a-425e-bb3e-386c8db8cd27
# ╟─0b6ec47b-3610-453f-b0c3-a99d5fb83ab4
# ╟─c1bf8bf8-fb05-4731-8140-c0f6332496bd
# ╟─60559545-8f85-4027-8471-4ce3d7cc952a
# ╠═3c4b1028-bcf7-45b0-bba2-dab5a63724d0
# ╠═f83d3db6-39bb-4f3c-9870-790b37741646
# ╠═785a2042-bb12-4a2a-8117-c2cfb396d29c
# ╟─fb55b1e8-eb70-45e9-9139-38ba64a50205
# ╟─0e1f6708-53d8-4cca-acf3-e7b2cac1dd7a
# ╟─de7e3ebb-f9d2-4830-9bec-f034d965ab76
# ╟─d2543921-868a-47d5-94f8-a965161ef8b4
# ╠═02e9bc9e-9d86-4b1f-adef-38a67d6ba2d1
# ╟─b12592bf-ec54-435f-b795-183cf84bcb4c
# ╠═91d634c4-e4ec-49b8-a035-5df0740f9fb0
# ╠═474d1a20-91ec-4b50-9caf-1991de12fb9e
# ╠═8d98f83f-cdb6-4856-b9a2-0e8afe39702f
# ╠═a25e4ec0-e8f2-436e-9f3a-c0c7da6f7fe2
# ╠═b03aeea2-4591-4d2a-af22-441e9bb2b571
# ╠═f2ed8d60-28e8-4fd4-8c08-88cf0e26ac1d
# ╠═8d6bda99-27e1-4993-9af5-816910d7bb21
# ╠═05d7024b-8a67-48c2-af78-527e9f1f8718
# ╠═8acc524d-692e-47f4-8c4d-ecbe6a8f4226
# ╠═493fea80-ad5a-4fb1-add5-ad81e262ea9f
# ╠═7c302268-c603-4e46-adb1-ba047429c2e0
# ╠═3f1b1f9c-2ee0-4c8b-832f-d984669450ba
# ╠═58086727-ef4a-4456-bfdf-f820432f4cad
# ╠═e00fc7cb-05a7-4215-b7f6-5bad647ab090
# ╠═4c62a111-112c-4b8d-979d-dea84e35232a
# ╟─2d56d9c0-b23a-4399-a86e-77edbf25d8b5
# ╠═8cb9eca6-2ebe-4faa-b6c6-5ea267e64dc9
# ╠═1174095b-a339-47b2-9f21-3f5b2c52b4eb
# ╠═21e0d0a9-0fb4-4bd5-93a3-517481b753d7
# ╠═3d088149-daa9-4664-9a06-7685d4cbb00e
# ╠═343f3226-a6e1-4bf5-a816-f92cb74cf434
# ╠═e9f30692-375c-49fc-bf0d-f045951ae11f
# ╠═2a41bd71-3796-4a21-b7ce-da0b09d043de
# ╟─21d22e3c-a56a-4c03-a9aa-226acea40c72
# ╠═2fd9c2e4-fd39-44b7-8d88-bddaf46ee26b
# ╠═8bb918eb-cdb6-4a15-a344-88ec22df92b2
# ╟─592b8b00-91a0-4ab4-b26b-9fc358978e9a
# ╠═df6493dd-5730-4993-b45d-495e62f34636
# ╠═72ad7501-5a9c-4e06-91e9-a03baff07a83
# ╠═78f03589-85e8-402d-9474-4fe2e5b8caf4
# ╠═49f1594c-c271-4981-9753-da88be755b23
# ╠═a782e454-470f-4373-ad93-2b73c0a69b8b
# ╠═7aa5561f-682a-4c15-b6d3-4c3ef8572da1
# ╠═701ebc68-2809-4544-be85-c4fadae6c391
# ╠═c863a4fc-270b-4473-beaa-d34ba2d719de
# ╟─ebf32f98-3ba0-4650-83ca-4d7a5b310c05
# ╠═bbf59177-36a8-480c-af39-1998b37032d2
# ╟─b886c106-b325-464b-8e71-6dc46350e48c
# ╟─f7db71f7-9bac-4aad-a636-5c60ae4c259b
# ╟─fb17d03b-880c-4be6-bd50-f31bdd3e2fd4
# ╟─9dc77012-f845-4a5b-a059-fd13c1308b3d
# ╟─f57e8544-204a-47a0-93c9-24a2982a0aed
# ╟─8af9d40f-8d50-4376-82e4-01981a54893f
# ╟─92175df7-b060-4f6f-8acb-265bdbe6713a
# ╟─0d9dd1ea-9ef8-460c-8d95-20a46eb1ebb6
# ╟─f8e4b341-1613-4993-a193-437f2adf4a21
# ╟─cfce959c-c171-4923-b5f9-29666d9c20e4
# ╟─fe0b7b47-39bb-45f8-873e-9a4a07d7d1b3
# ╟─b6da427f-0c72-4de7-bdca-b85688ac9bed
# ╠═b94f222a-0c35-49d7-94f2-a2347efb6b2d
# ╟─e7541e93-8b7a-4c26-901e-d887931b8103
# ╟─e858ca22-5ddd-45aa-8815-c3acc3f9741e
# ╟─ea5aa2da-f64f-44a7-8c41-746c9e0145e7
# ╟─fca2fa7d-99fe-44bd-a34a-1b1552a3ae48
# ╟─72e90532-02da-4d31-88c3-7766393d9f92
# ╟─0b07ac3c-d61e-43ff-acff-f907019835d6
# ╟─473302a3-2541-49a1-8f21-906bc2c9d279
# ╟─6a5425e0-7fbf-4e62-b839-bb34f00bd026
# ╟─a46da6b4-c8c9-4c1f-9bb6-33899bc63336
# ╟─34397405-e939-464c-a384-5d6c8fb9ca83
# ╟─9b81f986-ea3c-49f5-be04-b12b392aacd6
# ╟─7b2e6266-6dfe-47b5-8f28-3b9d76f57a42
# ╟─c8621274-00c9-4f67-b285-c624d1346527
# ╟─7d63e3eb-6132-439a-bc68-05e73cd5b068
# ╟─d4e7987c-94bf-4194-8e16-88c2581832d1
# ╟─4c71c3d6-1f9e-4812-9f36-2e65a43c8207
# ╟─b2b81585-2144-4b37-b0a7-286c862bb812
# ╟─0e9b78b0-208f-4b29-96c3-4b4129416bdd
# ╟─9833f763-ece1-4739-a6d1-0f3182d32d5d
# ╟─e18232e8-31e6-4fbf-8c96-5c36c9e8602b
# ╟─9d89d84d-f107-4a8d-be18-55af8edb546e
# ╟─ab2e501d-1c2d-48e1-be40-942cb4e66bc8
# ╟─833476e3-f1f1-4042-87e0-4c21ef37f7a7
# ╟─04e04706-cf7c-428d-8e3e-70643cf7ed60
# ╠═98087935-894b-408c-b708-02f05d97a96c
# ╟─2c88f409-6763-41d8-b1cc-168df855728b
# ╟─38bb2810-3a30-4a39-99dd-83518c093be2
# ╠═c23c4063-947e-436a-910b-c43afdb9c244
# ╠═26caf850-c46e-453c-aad9-c7bdb71f858b
# ╟─b957ac39-458f-48b6-8683-229d09ecaea4
# ╠═d9877f46-c461-4615-bb0c-0c7541238899
# ╠═93465fa0-83bc-413b-b46e-7ed0d3d0d1a7
# ╠═4226fab5-719e-45a5-ad3d-af712dd743ce
# ╠═022f7788-a2d0-4489-af0b-71692f800fdf
# ╠═5a155b3a-c1b6-4edb-8c5d-14cd490f9ce8
# ╠═019c4b03-caf9-48eb-b75d-882582a19d94
# ╠═5140cedf-cdd4-43a3-8778-3f926a9f31c0
# ╠═c7690425-1a71-4fdf-9b41-e7f4e5f13d5b
# ╠═fe96f553-ec82-4997-a0ca-72f26f4aa5dd
# ╠═34bda3ea-4238-4670-9f6b-9207f63047cb
# ╠═72522f83-2b3c-4e40-a1de-c3508d50fdab
# ╠═ea26ca76-102d-4011-9f8b-0917c021cc9f
# ╠═cf17cae1-69ea-4050-b302-cfc5d30d6235
# ╠═6fb4b458-6192-4000-8efd-f6fb02488861
# ╠═31f212b1-dac3-4c77-a340-dd2837fc8807
# ╟─4441d217-5f2e-4052-8fda-e11672ae1813
# ╠═37521853-b548-44b1-b4ea-df73eb2682d2
# ╠═cd026a7a-39f7-429e-95ac-a3be255f1f29
# ╟─e0756a3b-998a-4f8d-9ab3-df50422ce14b
# ╠═44ae02ac-6f9a-4af5-8901-3428da03af22
# ╠═b05ae27f-13df-4ae6-8189-1c21d86ee2df
# ╠═e3b8c1b9-cdb6-4a06-bc39-180f4bc42318
