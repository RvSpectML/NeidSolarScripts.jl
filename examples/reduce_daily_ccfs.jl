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

# ╔═╡ 42b4e6fe-a897-11eb-2676-f7fd96f35a22
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
	ccf_path = "/mnt/data_simons/NEID/v0.7develop20210501/outputs3"
	daily_filename = "daily_ccfs.jld2"
	max_days_to_use = 365
	order_idx_to_use_for_snr = 60
	order_snr_threshold_max = 36779.37219973966
	min_obs_in_day = 5
end;

# ╔═╡ 15741d2d-5fc2-4dc0-baf1-3ba3a72b5399
begin  # Parameters for picking which days to calculate RVs for
	order_snr_threshold = order_snr_threshold_max*0.94
	solar_hour_angle_threshold = 2.0
end;

# ╔═╡ d6219a57-7c1d-4b23-9b80-ae96f005b762
begin   # Parameters for building CCF template
	solar_hour_angle_threshold_template = 1.0
	order_snr_threshold_template = order_snr_threshold_max*0.96
end;

# ╔═╡ 60559545-8f85-4027-8471-4ce3d7cc952a
md"""
### Script
"""

# ╔═╡ 21d22e3c-a56a-4c03-a9aa-226acea40c72
md"""
# Results
## One Day
"""

# ╔═╡ fe0b7b47-39bb-45f8-873e-9a4a07d7d1b3
md"""
Order RV slope limits (lo, hi) $(@bind order_slope_lo Slider(-100.0:0.0; default=-20))  $(@bind order_slope_hi Slider(0.0:100.0; default=20.0))

"""

# ╔═╡ b6da427f-0c72-4de7-bdca-b85688ac9bed
md"""
## Combined RVs
#### Order indices to include in RV calculation
1 $(@bind use_ord_1 CheckBox(default=true))
2 $(@bind use_ord_2 CheckBox(default=false))
3 $(@bind use_ord_3 CheckBox(default=false))
4 $(@bind use_ord_4 CheckBox(default=true))
5 $(@bind use_ord_5 CheckBox(default=true))
6 $(@bind use_ord_6 CheckBox(default=true))
7 $(@bind use_ord_7 CheckBox(default=false))
8 $(@bind use_ord_8 CheckBox(default=true))
9 $(@bind use_ord_9 CheckBox(default=true))
10 $(@bind use_ord_10 CheckBox(default=true))
11 $(@bind use_ord_11 CheckBox(default=false))
12 $(@bind use_ord_12 CheckBox(default=true))
13 $(@bind use_ord_13 CheckBox(default=false))
14 $(@bind use_ord_14 CheckBox(default=false))
15 $(@bind use_ord_15 CheckBox(default=false))
16 $(@bind use_ord_16 CheckBox(default=false))
17 $(@bind use_ord_17 CheckBox(default=false))
18 $(@bind use_ord_18 CheckBox(default=true))
19 $(@bind use_ord_19 CheckBox(default=true))
20 $(@bind use_ord_20 CheckBox(default=true))

21 $(@bind use_ord_21 CheckBox(default=true))
22 $(@bind use_ord_22 CheckBox(default=true))
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
order_idx_for_rvs_manual = vcat(1,4:6,10:11,19:23);

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
$(@bind order_rv_rms_hi Slider(0.5:0.2:200.0; default=10))
"""

# ╔═╡ 9833f763-ece1-4739-a6d1-0f3182d32d5d
md"### Checking for Chromatic Effects"

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
datafiles = make_dataframe_of_ccf_files(ccf_path,daily_filename);

# ╔═╡ 26caf850-c46e-453c-aad9-c7bdb71f858b
function select_days_to_analyze_add_manifest_new(df::DataFrame; solar_hour_angle_threshold::Real = 3.0, order_snr_threshold::Real = order_snr_threshold_max*0.9, order_idx_to_use_for_snr::Integer = 60, min_obs_in_day::Integer = 5, max_days_to_use::Integer = 10)
	df[!,:nobs_in_file] .= 0
	#df[!,:nobs_use] .= 0
	df[!,:analyze_day] .= false
	df[!,:manifest] = fill(DataFrame(),size(df,1))
	
	num_days_to_use_found = 0
	for (i,day) in enumerate(eachrow(df))
		solar_info_fn = joinpath(dirname(day.filename),"..",day.datestr * "_solar_info.csv")
		println("Looking for ", solar_info_fn)
		@assert isfile(solar_info_fn) && filesize(solar_info_fn)>0
		df_solarinfo = CSV.read(solar_info_fn, DataFrame)
		manifest_all = load(day.filename,"manifest")
		df_combo = innerjoin(manifest_all, df_solarinfo, on = :bjd, makeunique = true)
		println("# Sizes: manifest= ", size(manifest_all,1), " solarinfo= ", size(df_solarinfo), " df_combo= ", size(df_combo,1))
		#orders_to_use = load(day.filename,"orders_to_use")
		#=
		manifest_use = manifest_all |> 
			@filter(abs(_.solar_hour_angle)<= solar_hour_angle_threshold) |> 
			@filter(_.order_snrs[order_idx_to_use_for_snr] >= order_snr_threshold ) |>
				DataFrame
		=#
		manifest_use = df_combo |> 
			@filter(abs(_.solar_hour_angle)<= solar_hour_angle_threshold) |> 
			#@filter(_.order_snrs[order_idx_to_use_for_snr] >= order_snr_threshold ) |>
			@filter( 0.995 <= _.pyrheliometer_flux_over_model <= 1.015 ) |>
			@filter( _.pyrheliometer_rms_flux <= 0.005 ) |>
				DataFrame
		
		day.nobs_in_file = size(manifest_all,1)
		#day.
		nobs_use = size(manifest_use,1)
		#day.manifest = manifest_all
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
#datafiles_use = select_days_to_analyze_add_manifest(datafiles,
datafiles_use = select_days_to_analyze_add_manifest_new(datafiles,
					max_days_to_use=max_days_to_use,
					solar_hour_angle_threshold=solar_hour_angle_threshold,
					order_snr_threshold=order_snr_threshold);

# ╔═╡ 433c64f5-3f0c-4442-a5ea-5d57e14a2038
size(datafiles_use)

# ╔═╡ 4ae8bce1-5890-4647-94e5-e0dfd14e9eb5
names(datafiles_use)

# ╔═╡ 188cf9ab-f839-490b-99fa-04611c355bfe
datafiles_use

# ╔═╡ 0e1f6708-53d8-4cca-acf3-e7b2cac1dd7a
begin   # Get size of order ccfs
	v_grid =  load(first(datafiles_use.filename),"v_grid")
	num_v = length(v_grid)
	ex_order_ccfs = load(first(datafiles_use.filename),"order_ccfs")
	num_order_idx = size(ex_order_ccfs,2)
end;

# ╔═╡ ebf32f98-3ba0-4650-83ca-4d7a5b310c05
md"""
## Orders RVs
Order index (min, max, step) 
$(@bind order_idx_lo NumberField(1:num_order_idx; default=4))  
$(@bind order_idx_hi NumberField(1:num_order_idx; default=14))
$(@bind order_idx_step NumberField(1:10; default=1))
Use order range? $(@bind use_order_range CheckBox(default=true))
Plot bad days? $(@bind plt_all_days CheckBox(default=false))
"""

# ╔═╡ c487c08a-93a1-42dc-b465-3268be5954df
datafiles_use[79,:manifest]

# ╔═╡ bcbd90c4-60b1-4c90-aa4b-69b0634e2fcd
datafiles_use.manifest[13]

# ╔═╡ 0b6ec47b-3610-453f-b0c3-a99d5fb83ab4
histogram(datafiles_use.manifest[1][!,:pyrheliometer_rms_flux],bins=range(0,stop=0.006,length=40))

# ╔═╡ c1bf8bf8-fb05-4731-8140-c0f6332496bd
histogram(datafiles_use.manifest[1][!,:pyrheliometer_flux_over_model],bins=range(0.995,stop=1.015,length=40))

# ╔═╡ b957ac39-458f-48b6-8683-229d09ecaea4
function select_days_to_analyze_add_manifest(df::DataFrame; solar_hour_angle_threshold::Real = 3.0, order_snr_threshold::Real = order_snr_threshold_max*0.9, order_idx_to_use_for_snr::Integer = 60, min_obs_in_day::Integer = 5, max_days_to_use::Integer = 10)
	df[!,:nobs_in_file] .= 0
	#df[!,:nobs_use] .= 0
	df[!,:analyze_day] .= false
	df[!,:manifest] = fill(DataFrame(),size(df,1))
	
	num_days_to_use_found = 0
	for (i,day) in enumerate(eachrow(df))
		manifest_all = load(day.filename,"manifest")
		#orders_to_use = load(day.filename,"orders_to_use")
		manifest_use = manifest_all |> 
			@filter(abs(_.solar_hour_angle)<= solar_hour_angle_threshold) |> 
			@filter(_.order_snrs[order_idx_to_use_for_snr] >= order_snr_threshold ) |>
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
function compute_mean_order_ccf(df::DataFrame; solar_hour_angle_threshold::Real = 1.0, order_snr_threshold::Real = order_snr_threshold_max*0.95, order_idx_to_use_for_snr::Integer = 60 )
	ex_order_ccfs = load(first(df.filename),"order_ccfs")
	mean_order_ccf = zeros(size(ex_order_ccfs,1),size(ex_order_ccfs,2))
	weights_order_ccf = zeros(size(ex_order_ccfs,1),size(ex_order_ccfs,2))
	for row in eachrow(df)
		sha_mask = abs.(row.manifest.solar_hour_angle).<= solar_hour_angle_threshold
	    snr_mask = map(i->row.manifest.order_snrs[i][order_idx_to_use_for_snr],1:length(row.manifest.order_snrs)) .>= order_snr_threshold
		mask = sha_mask .& snr_mask
		nobs_use = sum(mask)
		v_grid_this =  load(row.filename,"v_grid")
		@assert all(v_grid .== v_grid_this)
		order_ccfs = load(row.filename,"order_ccfs")
		order_ccf_vars = load(row.filename,"order_ccf_vars")
		@assert size(order_ccfs,1) == size(mean_order_ccf,1)
		@assert size(order_ccfs,2) == size(mean_order_ccf,2)
		num_obs = size(order_ccfs,3)
		weights = zeros(num_obs)
		for i in 1:num_v
			for j in 1:num_order_idx
				weights .= 1.0 ./ order_ccf_vars[i,j,:].^2
				mean_order_ccf[i,j] += sum(order_ccfs[i,j,mask] .* weights[mask])
				weights_order_ccf[i,j] += sum(weights[mask])
			end
		end
	end
	mean_order_ccf ./= weights_order_ccf
	return mean_order_ccf
end

# ╔═╡ de7e3ebb-f9d2-4830-9bec-f034d965ab76
mean_order_ccf = compute_mean_order_ccf(datafiles_use, solar_hour_angle_threshold=solar_hour_angle_threshold_template, order_snr_threshold=order_snr_threshold_template, order_idx_to_use_for_snr=order_idx_to_use_for_snr);

# ╔═╡ b12592bf-ec54-435f-b795-183cf84bcb4c
begin 
	mrv = []
	for ord_idx in 1:num_order_idx
		templ = vec(mean_order_ccf[:,ord_idx])
		if all(isnan.(templ)) 
			push!(mrv,missing) 
			continue
		end
		mrv_this_order = RVFromCCF.MeasureRvFromCCFTemplate(v_grid=v_grid,template=templ, frac_of_width_to_fit = 1.0, measure_width_at_frac_depth=0.5)
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

# ╔═╡ 022f7788-a2d0-4489-af0b-71692f800fdf
function find_good_day_idx(df::DataFrame; solar_hour_angle_threshold::Real = 1.0, order_snr_threshold::Real = order_snr_threshold_max*0.96, order_idx_to_use_for_snr::Integer = 60)
	argmax(
		df |> @map( sum((abs.(_.manifest.solar_hour_angle).<solar_hour_angle_threshold) .& map(i->_.manifest.order_snrs[i][order_idx_to_use_for_snr] > order_snr_threshold,1:size(_.manifest,1))) ) |> collect )
end

# ╔═╡ 8cb9eca6-2ebe-4faa-b6c6-5ea267e64dc9
obs_idx_good_day = find_good_day_idx(datafiles_use)

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
function summarize_day_orders(day_to_plt::Integer; df::DataFrame, orders, weights)
	local δrvs = df.manifest[day_to_plt][!,:Δv_diff_ext]
	local daily_means = map(ord->NaNMath.mean(df.rvs[day_to_plt][ord,:].+δrvs),orders)
	local daily_rmss = map(ord->NaNMath.std(df.rvs[day_to_plt][ord,:].+δrvs),orders)	
	local times = df.manifest[day_to_plt].bjd
	
	local daily_rmss_binned =  map(ord-> 
		NaNMath.std(bin_times_and_rvs_max_Δt(times=times, 
			rvs=df.rvs[day_to_plt][ord,:].+δrvs, 
				Δt_threshold= 6.0/(24*60))[2] ),orders)
	
	local shas = df.manifest[day_to_plt].solar_hour_angle ./24
	local tsha0 = roots(Polynomials.fit(times.-2.459e6,shas,1))[1].+2.459e6
	local daily_fits = map(ord->Polynomials.fit(times.-tsha0,df.rvs[day_to_plt][ord,:].+δrvs,1),orders)
	local daily_rv_sha0s = map(i->daily_fits[i][0],1:length(orders))
	local daily_slopes = map(i->daily_fits[i][1],1:length(orders))
	
	return Dict(:mean_order_rv=>daily_means,
		:rms_order_rv=>daily_rmss, :t_sha0=>tsha0, :Δt=>times.-tsha0, 
		:order_rv_slope=>daily_slopes, :order_rv_sha0=>daily_rv_sha0s,
		:rms_order_rv_binned=>daily_rmss_binned )
end

# ╔═╡ 4c62a111-112c-4b8d-979d-dea84e35232a
day_orders = DataFrame(summarize_day_orders.(1:size(datafiles_use,1),df=df_rvs,orders=1:num_order_idx, weights= order_weights))

# ╔═╡ 2d56d9c0-b23a-4399-a86e-77edbf25d8b5
median_order_rv = map(ord->NaNMath.median(map(i->day_orders.mean_order_rv[i][ord],(1:size(day_orders,1)))),1:num_order_idx)

# ╔═╡ 701ebc68-2809-4544-be85-c4fadae6c391
begin
	daily_rms_across_orders = map(obsid->NaNMath.std((day_orders.order_rv_sha0[obsid].-median_order_rv)[vcat(1,4:6,8:10,12)]),1:size(day_orders,1))
end

# ╔═╡ 6fb4b458-6192-4000-8efd-f6fb02488861
function summarize_day_total(day_to_plt::Integer; df::DataFrame, orders, weights)
	rv_tot = vec(sum(df.rvs[day_to_plt][orders,:].*weights[orders],dims=1) ./ 					   sum(weights[orders],dims=1))
	rv_tot .+= df.manifest[day_to_plt][!,:Δv_diff_ext]
	times = df.manifest[day_to_plt].bjd
	shas = df.manifest[day_to_plt].solar_hour_angle ./24
	tsha0 = roots(Polynomials.fit(times.-2.459e6,shas,1))[1].+2.459e6

	sha_mask = abs.(df.manifest[day_to_plt].solar_hour_angle).<= solar_hour_angle_threshold
	snr_mask = map(i->df.manifest[day_to_plt].order_snrs[i][order_idx_to_use_for_snr],1:length(df.manifest[day_to_plt].order_snrs)) .>= order_snr_threshold
	pyrh_mask = (0.995 .<= df.manifest[day_to_plt].pyrheliometer_flux_over_model .<= 1.015 ) .& ( df.manifest[day_to_plt].pyrheliometer_rms_flux .<= 0.006 )	
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

	Dict( :rv=>rv_tot, :rms_rv=>rms_rv,
		  :t_sha0=>tsha0, :Δt => times.-tsha0, 
		  :rv_slope => fit_slope, :rv_sha0 =>fit_constant,
		  :Δt_binned => t_binned.-tsha0, :rv_binned=>rvs_binned, :num_obs_in_bin=>num_obs_in_bin,
		  :rms_rv_binned=>rms_rv_binned, :rms_rv_binned_fit =>rms_rv_binned_fit,
		  :mask=>mask
		)
end

# ╔═╡ 343f3226-a6e1-4bf5-a816-f92cb74cf434
begin
	day_comb = DataFrame(summarize_day_total.(1:size(datafiles_use,1),df=df_rvs,orders=order_idx_for_rvs, weights= order_weights));
end;

# ╔═╡ 592b8b00-91a0-4ab4-b26b-9fc358978e9a
md"""
Day to plot
$(@bind obs_idx_plt NumberField(1:size(day_comb,1); default=obs_idx_good_day))  
"""

# ╔═╡ c863a4fc-270b-4473-beaa-d34ba2d719de
begin
	mask_decent_days = (length.(day_comb.rv_binned).>=20) #.&
		#(daily_rms_across_orders .< 6) #.&
	 	#(map(obsidx->median(map(ordidx->day_orders[obsidx,:rms_order_rv][ordidx],order_idx_for_rvs)),1:size(day_orders,1)) .< 6)
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

# ╔═╡ a79fb3f3-a07c-43b8-9be6-72e7347ff40e
size(days_to_plt)

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
	local orders_to_plt = use_order_range ? (order_idx_lo:order_idx_step:order_idx_hi) : order_idx_for_rvs
	local pal = palette(:berlin, length(orders_to_plt))
	for (i,ordidx) in enumerate(orders_to_plt)
		scatter!(plt,Date.(Dates.julian2datetime.(day_orders.t_sha0[days_to_plt])),	map(obsidx->day_orders[obsidx,:rms_order_rv][ordidx],1:size(day_orders,1)),ms=1.0,label=:none,c=pal[i])
		scatter!(plt,Date.(Dates.julian2datetime.(day_orders.t_sha0[days_to_plt])),	map(obsidx->day_orders[obsidx,:rms_order_rv][ordidx],days_to_plt),label=:none,c=pal[i])
	end
	xlabel!(plt,"Time (d)")
	ylabel!(plt,"RMS of Orders RVs within each day")
	ylims!(order_daily_rms_lo,order_daily_rms_hi)
	title!(plt,"Daily RMS of Order RVs")	
end

# ╔═╡ 9dc77012-f845-4a5b-a059-fd13c1308b3d
begin
	histogram(daily_rms_across_orders, 
		bins=range(order_daily_rms_lo,stop=order_daily_rms_hi,length=20),alpha=0.25, label="All usable days")
	histogram!(daily_rms_across_orders[days_to_plt], 
		bins=range(order_daily_rms_lo,stop=order_daily_rms_hi,length=20),alpha=0.5, label="Good days")
	xlabel!("RMS (m/s)")
	ylabel!("Count")
	title!("RMS in daily RV over orders used for RVs")
	xlims!(order_daily_rms_lo,order_daily_rms_hi)
end

# ╔═╡ f8e4b341-1613-4993-a193-437f2adf4a21
md"""
Order RV limits (lo, hi) $(@bind order_rv_lo Slider(-500.0:0.0; default=default_order_rv_lo))  $(@bind order_rv_hi Slider(0.0:500.0; default=default_order_rv_hi))
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
	local plt = plot()
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
	xlabel!("RMS (m/s)")
	ylabel!("Count")
	xlims!(order_rv_rms_lo,order_rv_rms_hi)
	plt
end

# ╔═╡ 5bc2c3e1-62dc-4f64-9339-87c6cbca8a8b
size(day_comb[44,:rv_binned])

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
	histogram(day_comb.rms_rv[days_to_plt], bins=daily_rv_rms_lo:(round((daily_rv_rms_hi-daily_rv_rms_lo)/40,digits=3)):daily_rv_rms_hi, alpha=0.2,label="Unbinned")
	histogram!(day_comb.rms_rv_binned[days_to_plt], bins=daily_rv_rms_lo:(round((daily_rv_rms_hi-daily_rv_rms_lo)/40,digits=3)):daily_rv_rms_hi,alpha=0.6, label="Binned")
	histogram!(day_comb.rms_rv_binned_fit[days_to_plt], bins=daily_rv_rms_lo:(round((daily_rv_rms_hi-daily_rv_rms_lo)/40,digits=3)):daily_rv_rms_hi, alpha=0.4, label="Binned - linear fit")
	
	xlims!(daily_rv_rms_lo,daily_rv_rms_hi)
	xlabel!("RMS (m/s)")
	ylabel!("Count")
end

# ╔═╡ 7b2e6266-6dfe-47b5-8f28-3b9d76f57a42
begin
	plot(legend=:none)
	scatter!(Date.(Dates.julian2datetime.(day_comb.t_sha0)),day_comb.rms_rv,ms=1.5,color=:blue)
	scatter!(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt])),day_comb.rms_rv[days_to_plt], color=:blue)
	scatter!(Date.(Dates.julian2datetime.(day_comb.t_sha0)),day_comb.rms_rv_binned, ms=2,color=:red)
	scatter!(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt])),day_comb.rms_rv_binned[days_to_plt], color=:red)
	title!("Daily RV RMS")
	xlabel!("Day")
	ylabel!("RMS RV within day (m/s)")
	xlims!(extrema(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt]))))
	ylims!(daily_rv_rms_lo,daily_rv_rms_hi)
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
	scatter!(plt,Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt])), Δrv[days_to_plt])
	local pred = long_term_daily_rv_fit.(day_comb.t_sha0)
	plot!(plt,Date.(Dates.julian2datetime.(day_comb.t_sha0)),pred)
	rms_resid = std(Δrv[days_to_plt][fit_mask].-pred[days_to_plt][fit_mask])
	rms_resid_0 = std(Δrv[days_to_plt][fit_mask])
	title!(plt,"RV fit at midday.  Linear slope: " * string(round(long_term_daily_rv_fit[1],digits=4)) * " m/s/d.\nRMS (fit): " * string(round(rms_resid,digits=3)) * "m/s"  * "  RMS (const): " * string(round(rms_resid_0,digits=3)) * "m/s")
	xlabel!(plt,"Day")
	ylabel!(plt,"RV at solar noon (m/s)")
	xlims!((minimum(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt].-1))),maximum(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt].+1)))))
	ylims!(plt,daily_rv_lo,daily_rv_hi)
	plt
end

# ╔═╡ 1c42a7ce-3e97-48b3-8ddd-89cab6177bfe
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
	histogram(day_comb.rv_slope[days_to_plt],bins=range(daily_slope_lo,stop=daily_slope_hi,length=30),alpha=0.5,label=:none)
	#local idx_good = length.(day_comb.rv_binned[days_to_plt]).>=25
	#histogram!(day_comb.rv_slope[days_to_plt][idx_good],bins=range(daily_slope_lo,stop=daily_slope_hi,length=30),alpha=0.5,label="Many points")
	title!("Distribution of Daily RV slopes")
	xlabel!("RV Slope (m/s/day)")
	xlims!(daily_slope_lo,daily_slope_hi)
end

# ╔═╡ e9f30692-375c-49fc-bf0d-f045951ae11f
begin
	day_farred = DataFrame(summarize_day_total.(1:size(datafiles_use,1),df=df_rvs,orders=vcat(36:40,46:50), weights= order_weights))
	day_red = DataFrame(summarize_day_total.(1:size(datafiles_use,1),df=df_rvs,orders=vcat(14:22,24:25), weights= order_weights))	
	day_blue = DataFrame(summarize_day_total.(1:size(datafiles_use,1),df=df_rvs,orders=vcat(1,4:6,8:10,12), weights= order_weights))
end;

# ╔═╡ 9d89d84d-f107-4a8d-be18-55af8edb546e
begin
	local plt = plot()
	local pfit = Polynomials.fit(day_comb.t_sha0[days_to_plt],Δrv[days_to_plt],1)
	Δrv_blue = day_blue.rv_sha0.-NaNMath.mean(day_blue.rv_sha0[days_to_plt])	
	Δrv_red =  day_red.rv_sha0.-NaNMath.mean(day_red.rv_sha0[days_to_plt])	
	Δrv_farred =  day_farred.rv_sha0.-NaNMath.mean(day_farred.rv_sha0[days_to_plt])	

	#scatter!(plt,Date.(Dates.julian2datetime.(day_blue.t_sha0[days_to_plt])), Δrv_blue[days_to_plt],color=:blue)
	scatter!(plt,Date.(Dates.julian2datetime.(day_red.t_sha0[days_to_plt])), Δrv_red[days_to_plt].-Δrv_blue[days_to_plt],color=:green, label="Red-Blue")
	scatter!(plt,Date.(Dates.julian2datetime.(day_farred.t_sha0[days_to_plt])), Δrv_farred[days_to_plt].-Δrv_blue[days_to_plt],color=:red, label="FarRed-Red")
	#rms_resid = std(Δrv[days_to_plt].-pfit.(day_comb.t_sha0[days_to_plt]))
	title!(plt,"ΔRV fit at midday") # \nLinear slope: " * string(round(pfit[1],digits=4)) * " m/s/d.  RMS about fit: " * string(round(rms_resid,digits=3)) * "m/s")
	xlabel!(plt,"Day")
	ylabel!(plt,"ΔRV at solar noon (m/s)")
	xlims!(extrema(Date.(Dates.julian2datetime.(day_comb.t_sha0[days_to_plt]))))
	ylims!(plt,-6,6)#daily_rv_lo,daily_rv_hi)
	plt
end

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
	ylims!(daily_slope_lo-2,daily_slope_hi+2)
end

# ╔═╡ 31f212b1-dac3-4c77-a340-dd2837fc8807
function plot_rvs_day(obs_idx::Integer)
	local plt = plot()
	sha_mask = abs.(datafiles_use.manifest[obs_idx].solar_hour_angle).<= solar_hour_angle_threshold
	 snr_mask = map(i->datafiles_use.manifest[obs_idx].order_snrs[i][order_idx_to_use_for_snr],1:length(datafiles_use.manifest[obs_idx].order_snrs)) .>= order_snr_threshold
	mask = sha_mask .& snr_mask
	plot!(plt,day_comb.Δt[obs_idx][mask].*24,day_comb.Δt[obs_idx][mask].*day_comb.rv_slope[obs_idx].+day_comb.rv_sha0[obs_idx], label=:none)
	scatter!(plt,day_comb.Δt[obs_idx][mask].*24,	day_comb.rv[obs_idx][mask], ms=1.5,label="Unbinned")
	(t_binned, rvs_binned) = bin_times_and_rvs_max_Δt(times=day_comb.Δt[obs_idx][mask],rvs=day_comb.rv[obs_idx][mask],Δt_threshold= 5/(24*60))
	scatter!(plt,t_binned.*24,	rvs_binned, label="Binned 5min")
	xlabel!(plt,"Δt from solar noon (hr)")
	ylabel!(plt,"RV (m/s)")
	title!(plt,datafiles_use.datestr[obs_idx] * ": RMS= " * string(round(day_comb.rms_rv_binned[obs_idx],digits=3)) * "m/s Slope= " * string(round(day_comb.rv_slope[obs_idx],digits=3)) * "m/s/d" )
	return plt
end

# ╔═╡ 8bb918eb-cdb6-4a15-a344-88ec22df92b2
plot_rvs_day(obs_idx_plt)

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
    max-width: 1030px;
    align-self: flex-start;
    margin-left: 20px;
}
"""

# ╔═╡ Cell order:
# ╟─87597e02-19c7-418a-8eb3-ad25695bb5bf
# ╟─ba3a32ee-9664-4279-b27b-0fed992e4d1c
# ╟─f64f9468-25ee-4799-b776-b1cfbe5b7f57
# ╠═c259acbf-da3d-4b73-af2b-8fc22f2548d7
# ╠═15741d2d-5fc2-4dc0-baf1-3ba3a72b5399
# ╠═d6219a57-7c1d-4b23-9b80-ae96f005b762
# ╟─60559545-8f85-4027-8471-4ce3d7cc952a
# ╠═3c4b1028-bcf7-45b0-bba2-dab5a63724d0
# ╠═785a2042-bb12-4a2a-8117-c2cfb396d29c
# ╠═433c64f5-3f0c-4442-a5ea-5d57e14a2038
# ╠═4ae8bce1-5890-4647-94e5-e0dfd14e9eb5
# ╠═188cf9ab-f839-490b-99fa-04611c355bfe
# ╟─0e1f6708-53d8-4cca-acf3-e7b2cac1dd7a
# ╟─de7e3ebb-f9d2-4830-9bec-f034d965ab76
# ╟─b12592bf-ec54-435f-b795-183cf84bcb4c
# ╠═91d634c4-e4ec-49b8-a035-5df0740f9fb0
# ╠═e00fc7cb-05a7-4215-b7f6-5bad647ab090
# ╟─8cb9eca6-2ebe-4faa-b6c6-5ea267e64dc9
# ╠═4c62a111-112c-4b8d-979d-dea84e35232a
# ╟─2d56d9c0-b23a-4399-a86e-77edbf25d8b5
# ╟─1174095b-a339-47b2-9f21-3f5b2c52b4eb
# ╠═343f3226-a6e1-4bf5-a816-f92cb74cf434
# ╠═e9f30692-375c-49fc-bf0d-f045951ae11f
# ╟─21d22e3c-a56a-4c03-a9aa-226acea40c72
# ╟─8bb918eb-cdb6-4a15-a344-88ec22df92b2
# ╟─592b8b00-91a0-4ab4-b26b-9fc358978e9a
# ╟─701ebc68-2809-4544-be85-c4fadae6c391
# ╠═c863a4fc-270b-4473-beaa-d34ba2d719de
# ╠═5bc2c3e1-62dc-4f64-9339-87c6cbca8a8b
# ╠═a79fb3f3-a07c-43b8-9be6-72e7347ff40e
# ╟─ebf32f98-3ba0-4650-83ca-4d7a5b310c05
# ╟─fb17d03b-880c-4be6-bd50-f31bdd3e2fd4
# ╟─bbf59177-36a8-480c-af39-1998b37032d2
# ╟─f7db71f7-9bac-4aad-a636-5c60ae4c259b
# ╟─9dc77012-f845-4a5b-a059-fd13c1308b3d
# ╟─8af9d40f-8d50-4376-82e4-01981a54893f
# ╟─92175df7-b060-4f6f-8acb-265bdbe6713a
# ╟─0d9dd1ea-9ef8-460c-8d95-20a46eb1ebb6
# ╟─f8e4b341-1613-4993-a193-437f2adf4a21
# ╟─cfce959c-c171-4923-b5f9-29666d9c20e4
# ╟─fe0b7b47-39bb-45f8-873e-9a4a07d7d1b3
# ╟─b6da427f-0c72-4de7-bdca-b85688ac9bed
# ╠═b94f222a-0c35-49d7-94f2-a2347efb6b2d
# ╟─e7541e93-8b7a-4c26-901e-d887931b8103
# ╠═e858ca22-5ddd-45aa-8815-c3acc3f9741e
# ╟─1c42a7ce-3e97-48b3-8ddd-89cab6177bfe
# ╟─0b07ac3c-d61e-43ff-acff-f907019835d6
# ╟─473302a3-2541-49a1-8f21-906bc2c9d279
# ╟─34397405-e939-464c-a384-5d6c8fb9ca83
# ╟─9b81f986-ea3c-49f5-be04-b12b392aacd6
# ╟─7b2e6266-6dfe-47b5-8f28-3b9d76f57a42
# ╟─c8621274-00c9-4f67-b285-c624d1346527
# ╟─7d63e3eb-6132-439a-bc68-05e73cd5b068
# ╟─fca2fa7d-99fe-44bd-a34a-1b1552a3ae48
# ╟─72e90532-02da-4d31-88c3-7766393d9f92
# ╠═c487c08a-93a1-42dc-b465-3268be5954df
# ╟─d4e7987c-94bf-4194-8e16-88c2581832d1
# ╟─4c71c3d6-1f9e-4812-9f36-2e65a43c8207
# ╟─b2b81585-2144-4b37-b0a7-286c862bb812
# ╟─0e9b78b0-208f-4b29-96c3-4b4129416bdd
# ╟─9833f763-ece1-4739-a6d1-0f3182d32d5d
# ╟─9d89d84d-f107-4a8d-be18-55af8edb546e
# ╟─833476e3-f1f1-4042-87e0-4c21ef37f7a7
# ╟─2c88f409-6763-41d8-b1cc-168df855728b
# ╟─38bb2810-3a30-4a39-99dd-83518c093be2
# ╠═c23c4063-947e-436a-910b-c43afdb9c244
# ╠═26caf850-c46e-453c-aad9-c7bdb71f858b
# ╠═bcbd90c4-60b1-4c90-aa4b-69b0634e2fcd
# ╠═0b6ec47b-3610-453f-b0c3-a99d5fb83ab4
# ╠═c1bf8bf8-fb05-4731-8140-c0f6332496bd
# ╠═b957ac39-458f-48b6-8683-229d09ecaea4
# ╟─d9877f46-c461-4615-bb0c-0c7541238899
# ╟─93465fa0-83bc-413b-b46e-7ed0d3d0d1a7
# ╟─022f7788-a2d0-4489-af0b-71692f800fdf
# ╟─34bda3ea-4238-4670-9f6b-9207f63047cb
# ╠═72522f83-2b3c-4e40-a1de-c3508d50fdab
# ╠═6fb4b458-6192-4000-8efd-f6fb02488861
# ╟─31f212b1-dac3-4c77-a340-dd2837fc8807
# ╟─4441d217-5f2e-4052-8fda-e11672ae1813
# ╟─42b4e6fe-a897-11eb-2676-f7fd96f35a22
# ╟─e0756a3b-998a-4f8d-9ab3-df50422ce14b
# ╟─44ae02ac-6f9a-4af5-8901-3428da03af22
# ╠═b05ae27f-13df-4ae6-8189-1c21d86ee2df
