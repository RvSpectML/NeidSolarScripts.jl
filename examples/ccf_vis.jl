### A Pluto.jl notebook ###
# v0.14.4

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
	using JLD2, FileIO
	using PlutoUI
	using Plots, ColorSchemes
	using Statistics, NaNMath
	using Polynomials
end

# ╔═╡ 5313b674-a897-11eb-3820-218f26a14d4d
begin
	ccf_path = ".."
	ccf_fn = "daily_ccfs.jld2"
	@assert isfile(joinpath(ccf_path,ccf_fn))
	@assert filesize(joinpath(ccf_path,ccf_fn)) > 0
end

# ╔═╡ 37a6c623-e901-4f8b-8989-835958471549
begin
	data = load(joinpath(ccf_path,ccf_fn))
	num_obs = size(data["order_ccfs"],3)
	num_order_idx = size(data["order_ccfs"],2)
	num_v = size(data["order_ccfs"],1)
	num_pix = size(data["mean_clean_flux"],1)
	num_pix_continuum = size(data["mean_clean_flux_continuum_normalized"],1)
	date_str = match(r"neidL1_(\d+)T\d+\.fits",data["manifest"].Filename[1]).captures[1]
	size(data["order_ccfs"])
end

# ╔═╡ d1fe0408-cccd-45c9-a5ef-c13a42386e44
begin
	(vmin_default, vmax_default) = extrema(data["v_grid"])
	order_ccfs_norm = data["order_ccfs"]./maximum(data["order_ccfs"],dims=1)
	order_ccfs_mean = mean(order_ccfs_norm,dims=3)
	ccf_order_depths = 1.0 .- minimum(order_ccfs_norm,dims=1)./maximum(order_ccfs_norm,dims=1) 
	ccf_mean_order_depth = vec(mean(ccf_order_depths,dims=2))
end;

# ╔═╡ 609e3d45-e1a5-4060-97a8-5b668f8137c1
md"""
### Pixel ranges for figures plotting mean flux.
Colum  (min, max) $(@bind xmin NumberField(1:num_pix; default=1))  $(@bind xmax NumberField(1:num_pix; default=num_pix))
Show uncertainties?   
Raw $(@bind plot_raw_errs CheckBox(default=false))  
Continuum $(@bind plot_continum_errs CheckBox(default=false))  
"""

# ╔═╡ 28aced5a-3100-4523-9f1a-f7b14016791b
md"""
#### Order index & Observation for plotting single CCF.

Order index  $(@bind order_idx NumberField(1:num_order_idx; default=floor(Int64,median(1:num_order_idx))))
Obs index  $(@bind obs_idx NumberField(1:num_obs; default=1))
"""

# ╔═╡ abbc4880-fedd-4bdd-a1db-b94e8a4d5654
md" vmin $(@bind vmin_kms NumberField(-100:100; default=-30))
vmax $(@bind vmax_kms NumberField(-100:100; default=30))
Ready? $(@bind ready_to_plot CheckBox(default=false))  
"

# ╔═╡ 6c469209-c7e0-41f2-97e3-db587b894e56
begin
	vmin = vmin_kms *1000;
	vmax = vmax_kms *1000;
end

# ╔═╡ ab0fc462-a897-11eb-1727-e71788522c02
if !all(isnan.(data["order_ccfs"][:,order_idx,obs_idx]))
  plt1 = plot(data["v_grid"],data["order_ccfs"][:,order_idx,obs_idx], ms=1.2,legend=:none)
  xlims!(plt1,vmin,vmax)
	title!("CCF Order Idx " * string(data["orders_to_use"][order_idx]) * " Obs " * string(obs_idx))
  #ylims!(plt1,extrema(f_filtered[λmin.<=lambda.<=λmax]))
  plt1
end

# ╔═╡ c403d4ff-3035-4bfc-9b9e-3f33ed0a6b7a
md"""
### Ranges for figures plotting multiple CCFs.
Order index (min, max, step) $(@bind order_idx_lo NumberField(1:num_order_idx; default=1))  $(@bind order_idx_hi NumberField(1:num_order_idx; default=num_order_idx))
$(@bind order_idx_step NumberField(1:num_order_idx; default=1))
"""

# ╔═╡ d84af1b1-1023-4493-8b0f-90c9030be4bd
md"Scale CCF depths? $(@bind scale_order_ccf_depths CheckBox(default=false))  "

# ╔═╡ 7dc1bf77-174b-4614-8aea-15cad25c0335
begin
  plt2 = plot()
  pal2 = palette(:berlin, length(order_idx_lo:order_idx_step:order_idx_hi))
  for (i,ord_idx) in enumerate(order_idx_lo:order_idx_step:order_idx_hi)
		if all(isnan.(order_ccfs_norm[:,ord_idx,obs_idx])) continue end
  		y = scale_order_ccf_depths ? 1. .- (1.0 .-order_ccfs_norm[:,ord_idx,obs_idx])./ccf_order_depths[1,ord_idx,obs_idx] :
		order_ccfs_norm[:,ord_idx,obs_idx]
		#plot!(plt2,data["v_grid"],order_ccfs_norm[:,ord_idx,obs_idx], ms=1.2,color=pal2[i],legend=:none)
		plot!(plt2,data["v_grid"],y, ms=1.2,color=pal2[i],legend=:none)
  end
  xlims!(plt2,vmin,vmax)
  xlabel!(plt2,"v (m/s)")
  ylabel!(plt2,scale_order_ccf_depths ? "Scaled CCFs" : "Normalized CCFs")
  title!(plt2,"CCF Order Idxs " * string(data["orders_to_use"][order_idx_lo]) * "-" * string(data["orders_to_use"][order_idx_hi]) * " Obs " * string(obs_idx))
  #ylims!(plt2,extrema(f_filtered[λmin.<=lambda.<=λmax]))
  plt2
end

# ╔═╡ ceae0ff5-f650-4dbf-81e5-26864eb1adcc
begin
	local z = scale_order_ccf_depths ? 1. .- (1.0 .-order_ccfs_norm[:,order_idx_lo:order_idx_hi,obs_idx])./reshape(ccf_order_depths[1,order_idx_lo:order_idx_hi,obs_idx],1,length(order_idx_lo:order_idx_hi)) :
		order_ccfs_norm[:,order_idx_lo:order_idx_hi,obs_idx]
	plt5 = heatmap(data["v_grid"],data["orders_to_use"][order_idx_lo:order_idx_hi],z')
	xlims!(plt5,vmin,vmax)
end

# ╔═╡ df7571c9-8b3f-4503-8bbe-42a77f0f9c0e
md"""Obs index range (min, max, step) $(@bind obs_idx_lo NumberField(1:num_obs; default=1))    
$(@bind obs_idx_hi NumberField(1:num_obs; default=num_obs))
$(@bind obs_idx_step NumberField(1:num_obs; default=1))
"""

# ╔═╡ 468fe2c8-7c4c-45f2-9577-a2cecba64d2d
begin
  pal3 = palette(:berlin, length(obs_idx_lo:obs_idx_step:obs_idx_hi))
  plt3 = plot()
  for (i,obs_id) in enumerate(obs_idx_lo:obs_idx_step:obs_idx_hi)
		if all(isnan.(order_ccfs_norm[:,order_idx,i])) continue end
		plot!(plt3,data["v_grid"],order_ccfs_norm[:,order_idx,i].-order_ccfs_mean[:,order_idx], ms=1.2,c=pal3[i],legend=:none)
	end
  xlims!(plt3,vmin,vmax)
  #ylims!(plt2,extrema(f_filtered[λmin.<=lambda.<=λmax]))
  xlabel!(plt3,"v (m/s)")
  ylabel!(plt3,"Normalized CCF - Time ave")
  title!(plt3,"CCF Order Idx " * string(order_idx) * "   Obs " * string(obs_idx_lo) * "-" * string(obs_idx_hi) )
  plt3
end

# ╔═╡ 0e4e2de1-d706-4019-bece-e24367b38d33
begin
	plt4 = heatmap(data["v_grid"],obs_idx_lo:obs_idx_hi,(order_ccfs_norm[:,order_idx,obs_idx_lo:obs_idx_step:obs_idx_hi].-order_ccfs_mean[:,order_idx])')
	xlabel!(plt4,"v (m/s)")
	ylabel!(plt4,"Observation number")
end

# ╔═╡ ed5e02d9-cfba-4ee6-85f2-3678cb1899f9

begin
	if eltype(data["normalization_anchors"][1]) <:AbstractVector
	num_orders_continuum = length(first(data["normalization_anchors"]))
	anchor_density = zeros(num_pix,num_orders_continuum)
	sigma_pix = 3
	for i in 1:num_obs
		for j in 1:num_orders_continuum
			anchors_here = data["normalization_anchors"][i][j]
			n_anchors_here = length(anchors_here)
			for k in 1:n_anchors_here
				anchor_density[:,j] .+= map(l->exp(-0.5*(anchors_here[k]-l)^2/sigma_pix^2),1:num_pix)./(sigma_pix*sqrt(2π))
			end
		end
	end
	end
end;

# ╔═╡ d251eb17-a164-4a36-90f3-d1d5aa87405c
begin
	plt = plot(legend=:bottomleft)
	#xmin = 1800;  xmax = 9500;
	pixlo = searchsortedfirst(1:size(data["mean_clean_flux_continuum_normalized"],1),xmin)
	pixhi = searchsortedlast(1:size(data["mean_clean_flux_continuum_normalized"],1),xmax)
	yscale_raw = NaNMath.maximum(data["mean_clean_flux"][pixlo:pixhi,data["orders_to_use"][order_idx]])
	yscale_sed = NaNMath.maximum(data["mean_clean_flux_sed_normalized"][pixlo:pixhi,data["orders_to_use"][order_idx]])
	yscale_continuum = NaNMath.maximum(data["mean_clean_flux_continuum_normalized"][pixlo:pixhi,data["orders_to_use"][order_idx]])

	plot!(plt,pixlo:pixhi,
		data["mean_clean_flux_sed_normalized"][pixlo:pixhi,data["orders_to_use"][order_idx]]./yscale_sed, color=:green, label="SED Normalized")
	plot!(plt,pixlo:pixhi,
			data["mean_clean_flux_continuum_normalized"][pixlo:pixhi,data["orders_to_use"][order_idx]]./yscale_continuum, color=:red, label="Continuum Normalized")
	if plot_continum_errs 
		plot!(plt,pixlo:pixhi,
			data["mean_clean_flux_continuum_normalized"][pixlo:pixhi,data["orders_to_use"][order_idx]]./yscale_continuum.+sqrt.(data["mean_clean_var_continuum_normalized"][pixlo:pixhi,data["orders_to_use"][order_idx]])./yscale_continuum,color=:red, ls=:dot,label=:none)
		plot!(plt,pixlo:pixhi,
			data["mean_clean_flux_continuum_normalized"][pixlo:pixhi,data["orders_to_use"][order_idx]]./yscale_continuum.-sqrt.(data["mean_clean_var_continuum_normalized"][pixlo:pixhi,data["orders_to_use"][order_idx]])./yscale_continuum,color=:red, ls=:dot,label=:none)
	else
	end
	plot!(plt,pixlo:pixhi,data["mean_clean_flux"][pixlo:pixhi,data["orders_to_use"][order_idx]]./yscale_raw, color=:blue, label="Raw")
	if plot_raw_errs
		plot!(plt,pixlo:pixhi,data["mean_clean_flux"][pixlo:pixhi,data["orders_to_use"][order_idx]]./yscale_raw.+sqrt.(data["mean_clean_var"][pixlo:pixhi,data["orders_to_use"][order_idx]])./yscale_raw
			, color=:blue, ls=:dot, label=:none)
		plot!(plt,pixlo:pixhi,data["mean_clean_flux"][pixlo:pixhi,data["orders_to_use"][order_idx]]./yscale_raw.-sqrt.(data["mean_clean_var"][pixlo:pixhi,data["orders_to_use"][order_idx]])./yscale_raw
			, color=:blue, ls=:dot, label=:none)
	end
	if eltype(data["normalization_anchors"][1]) <:AbstractVector
		plot!(plt,anchor_density[:,order_idx]./maximum(anchor_density[:,order_idx]),color=:orange,label="Anchor density")
	else
		anchor_idx = data["normalization_anchors"][order_idx]
		if length(anchor_idx) >0
			scatter!(plt,anchor_idx,ones(length(anchor_idx)),color=:magenta,label="Anchors")
			scatter!(plt,anchor_idx,data["mean_clean_flux"][anchor_idx,data["orders_to_use"][order_idx]]./yscale_raw,color=:blue,ms=2,label=:none)
		end
	end
	xlims!(pixlo,pixhi)
	#ylims!(0.0,1.1)
	plt
end

# ╔═╡ eff49841-4793-486a-8ed0-1ebe6deb4759
begin 
	anchors_merged = []
	if eltype(data["normalization_anchors"][1]) <:AbstractVector
	threshold = maximum(anchor_density[:,order_idx])/2
	for i in 2:(num_pix-1)
		if anchor_density[i,order_idx] > threshold 
			if (anchor_density[i,order_idx]>anchor_density[i+1,order_idx]) && (anchor_density[i,order_idx]>anchor_density[i-1,order_idx])
				push!(anchors_merged, i)
			end
		end
	end
	end
	anchors_merged
end

# ╔═╡ 3e6f6f21-878d-4618-9865-c1e656bc3104
begin
	if eltype(data["normalization_anchors"][1]) <:AbstractVector
	plt12 = scatter(anchor_density[:,order_idx])
	xlabel!(plt12,"Pixel")
	ylabel!(plt12,"Anchor density")
	#xlims!(plt12,900,1100)
	plt12
	end
end

# ╔═╡ c7a9d874-d823-45ff-8214-da400d01f383
collect(keys(data))

# ╔═╡ e3a2d62e-2d0f-4306-bb26-474445ce894e
data["orders_to_use"]

# ╔═╡ 1d3b78eb-983a-4f74-b064-e88a387dc5a1
data["mean_clean_flux"]

# ╔═╡ e0d97999-78cc-4e89-9534-66ca44dffee5
data["mean_clean_flux_sed_normalized"]

# ╔═╡ 4e718a48-6b04-441d-a0f6-17af907780ed
data["mean_clean_flux_continuum_normalized"]

# ╔═╡ c15b050a-f0de-4079-8fcb-912630722788
all(isnan.(data["mean_clean_flux_continuum_normalized"][1,:]))

# ╔═╡ a0103542-f745-4913-a8cf-a3bddd1631e4
length(data["mean_clean_flux_continuum_normalized"][:,1])

# ╔═╡ d9782c42-a2f7-41ae-a324-0cc86a56f393
data["order_ccfs"]

# ╔═╡ 9bafebc2-daae-42ea-9f68-4020a590a030
if eltype(data["normalization_anchors"][1]) <:AbstractVector
	map(i->data["normalization_anchors"][i][21],1:length(data["normalization_anchors"]))
end

# ╔═╡ 2735c043-d57d-4782-8ad1-a1254e373fe0
keys(data)

# ╔═╡ 82f4ceb0-617a-434e-97de-c33d500958f4
begin
	mrv = map(ord_idx->RVFromCCF.MeasureRvFromCCFTemplate(v_grid=data["v_grid"],template=vec(mean(data["order_ccfs"],dims=3)[:,ord_idx,:]), frac_of_width_to_fit = 1.0, measure_width_at_frac_depth=0.5),1:num_order_idx)
	begin
		rv_matrix = zeros(num_order_idx, num_obs)
		σ_rv_matrix = zeros(num_order_idx, num_obs)
		for obs_idx in 1:num_obs
			for ord_idx in 1:num_order_idx
				(rv, σ_rv) = mrv[ord_idx](data["v_grid"],data["order_ccfs"][:,ord_idx,obs_idx],data["order_ccf_vars"][:,ord_idx,obs_idx]) 
				rv_matrix[ord_idx,obs_idx] = rv
				σ_rv_matrix[ord_idx,obs_idx] = σ_rv
			end
		end
	end
end

# ╔═╡ c7ec68a3-6a2b-4cad-b1ff-2cb9b23cc499
begin
	scatter(σ_rv_matrix[:,1])
	scatter!(map(i->std(rv_matrix[i,:]),1:21))
end

# ╔═╡ e7541e93-8b7a-4c26-901e-d887931b8103
order_idx_for_rvs = vcat(1:5,6,8,9,10,17,21)
#order_idx_for_rvs = vcat(1:11,19:21)
#order_idx_for_rvs = vcat(1:5,8:10,19:21)

# ╔═╡ 1a8e2e8d-97a9-4634-9121-0f33266227f8
scatter(data["orders_to_use"][1:21],map(i->std(rv_matrix[i,:]),1:21))

# ╔═╡ e7906830-5dfd-42ea-ad8e-948b7a3138db
begin
	rvs = map(i->NaNMath.sum(rv_matrix[order_idx_for_rvs,i]./σ_rv_matrix[order_idx_for_rvs,i].^2)./NaNMath.sum(1.0 ./σ_rv_matrix[order_idx_for_rvs,i].^2),1:num_obs)
	σ_rvs = sqrt.(map(i->1.0 ./NaNMath.sum(1.0 ./σ_rv_matrix[order_idx_for_rvs,i].^2),1:num_obs))
end;

# ╔═╡ 9d38c13c-fbba-4fa2-b528-feae03551f54
(t_binned, rvs_binned) = bin_times_and_rvs_max_Δt(times=data["manifest"].bjd,rvs=rvs,Δt_threshold= 5/(24*60))

# ╔═╡ aac89b26-bd5d-402e-9b1c-4ab6a70dfe51
(NaNMath.std(rvs), NaNMath.mean(σ_rvs))

# ╔═╡ 2f9d86f3-ce3f-4a9c-a624-e0967a0ada70
std(rvs_binned)

# ╔═╡ 92517b16-e8f2-4c16-847f-350c93aea60f


# ╔═╡ 21d940fc-c258-4f01-a88b-71ff2109f80d
daily_slope = Polynomials.fit(t_binned,rvs_binned,1)[1]

# ╔═╡ 493b5c40-214a-41cd-aaaf-02609704ad49
daily_fit = Polynomials.fit(data["manifest"].bjd.-mean(data["manifest"].bjd),rvs,1)

# ╔═╡ 2921145f-df6a-4b14-9faf-46a27e8f9057
#std(rvs.-daily_fit.(data["manifest"].bjd.-mean(data["manifest"].bjd))), 
binned_rms_minus_linear = std(rvs_binned.-daily_fit.(t_binned.-mean(data["manifest"].bjd)))

# ╔═╡ facd440e-0e06-4824-9da0-cc7dd549db91
begin
	plt11 = scatter(data["manifest"].bjd.-mean(data["manifest"].bjd),rvs,yerr=σ_rvs, ms=2.0,label=:none)
	scatter!(plt11,t_binned.-mean(data["manifest"].bjd),rvs_binned, ms=4,label="Binned")
	plot!(plt11,data["manifest"].bjd.-mean(data["manifest"].bjd),daily_fit.(data["manifest"].bjd.-mean(data["manifest"].bjd)),label="Daily slope = " * string(round(Polynomials.fit(t_binned,rvs_binned,1)[1]/24,digits=3)) * "m/s/hr")
	xlabel!(plt11,"Time (d)")
	ylabel!(plt11,"RV (m/s)")
	title!(plt11,date_str * "  RMS_5min-line = " * string(round(binned_rms_minus_linear,digits=3)) * "m/s")
	plt11
end

# ╔═╡ 8c9d4dbf-c36e-4eee-b49f-ed0be110b55b
scatter(t_binned,rvs_binned)

# ╔═╡ 01d28d16-6acc-4209-893f-cd8cb730922d
begin
	scatter(rv_matrix[1,:])
	scatter!(rv_matrix[2,:])
	scatter!(rv_matrix[3,:])
	scatter!(rv_matrix[4,:])
	scatter!(rv_matrix[5,:])
	scatter!(rv_matrix[6,:])
	scatter!(rv_matrix[8,:])
	scatter!(rv_matrix[9,:])
	scatter!(rv_matrix[10,:])
end

# ╔═╡ Cell order:
# ╠═42b4e6fe-a897-11eb-2676-f7fd96f35a22
# ╠═5313b674-a897-11eb-3820-218f26a14d4d
# ╠═37a6c623-e901-4f8b-8989-835958471549
# ╠═d1fe0408-cccd-45c9-a5ef-c13a42386e44
# ╟─609e3d45-e1a5-4060-97a8-5b668f8137c1
# ╟─d251eb17-a164-4a36-90f3-d1d5aa87405c
# ╟─28aced5a-3100-4523-9f1a-f7b14016791b
# ╟─abbc4880-fedd-4bdd-a1db-b94e8a4d5654
# ╟─6c469209-c7e0-41f2-97e3-db587b894e56
# ╟─ab0fc462-a897-11eb-1727-e71788522c02
# ╟─c403d4ff-3035-4bfc-9b9e-3f33ed0a6b7a
# ╟─d84af1b1-1023-4493-8b0f-90c9030be4bd
# ╟─7dc1bf77-174b-4614-8aea-15cad25c0335
# ╟─ceae0ff5-f650-4dbf-81e5-26864eb1adcc
# ╟─df7571c9-8b3f-4503-8bbe-42a77f0f9c0e
# ╟─468fe2c8-7c4c-45f2-9577-a2cecba64d2d
# ╟─0e4e2de1-d706-4019-bece-e24367b38d33
# ╟─ed5e02d9-cfba-4ee6-85f2-3678cb1899f9
# ╟─eff49841-4793-486a-8ed0-1ebe6deb4759
# ╟─3e6f6f21-878d-4618-9865-c1e656bc3104
# ╠═c7a9d874-d823-45ff-8214-da400d01f383
# ╠═e3a2d62e-2d0f-4306-bb26-474445ce894e
# ╠═1d3b78eb-983a-4f74-b064-e88a387dc5a1
# ╠═e0d97999-78cc-4e89-9534-66ca44dffee5
# ╠═4e718a48-6b04-441d-a0f6-17af907780ed
# ╠═c15b050a-f0de-4079-8fcb-912630722788
# ╠═a0103542-f745-4913-a8cf-a3bddd1631e4
# ╠═d9782c42-a2f7-41ae-a324-0cc86a56f393
# ╠═9bafebc2-daae-42ea-9f68-4020a590a030
# ╠═2735c043-d57d-4782-8ad1-a1254e373fe0
# ╟─82f4ceb0-617a-434e-97de-c33d500958f4
# ╠═c7ec68a3-6a2b-4cad-b1ff-2cb9b23cc499
# ╠═e7541e93-8b7a-4c26-901e-d887931b8103
# ╠═1a8e2e8d-97a9-4634-9121-0f33266227f8
# ╟─e7906830-5dfd-42ea-ad8e-948b7a3138db
# ╠═9d38c13c-fbba-4fa2-b528-feae03551f54
# ╠═aac89b26-bd5d-402e-9b1c-4ab6a70dfe51
# ╠═2f9d86f3-ce3f-4a9c-a624-e0967a0ada70
# ╠═92517b16-e8f2-4c16-847f-350c93aea60f
# ╠═21d940fc-c258-4f01-a88b-71ff2109f80d
# ╠═493b5c40-214a-41cd-aaaf-02609704ad49
# ╠═2921145f-df6a-4b14-9faf-46a27e8f9057
# ╠═facd440e-0e06-4824-9da0-cc7dd549db91
# ╠═8c9d4dbf-c36e-4eee-b49f-ed0be110b55b
# ╠═01d28d16-6acc-4209-893f-cd8cb730922d
