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
	using JLD2, FileIO
	using PlutoUI
	using Plots, ColorSchemes
	using Statistics
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
	size(data["order_ccfs"])
end

# ╔═╡ d1fe0408-cccd-45c9-a5ef-c13a42386e44
begin
	(vmin_default, vmax_default) = extrema(data["v_grid"])
	order_ccfs_norm = data["order_ccfs"]./maximum(data["order_ccfs"],dims=1)
	order_ccfs_mean = #reshape(
		mean(order_ccfs_norm,dims=3) #,num_v,num_order_idx)

	ccf_order_depths = #reshape(
		1.0 .- minimum(order_ccfs_norm,dims=1)./maximum(order_ccfs_norm,dims=1) #,num_order_idx,num_obs)
	ccf_mean_order_depth = vec(mean(ccf_order_depths,dims=2))
end;

# ╔═╡ 28aced5a-3100-4523-9f1a-f7b14016791b
md"""
#### Order index & Observation for plotting single CCF.

Order index  $(@bind order_idx NumberField(1:num_order_idx; default=1))
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
	title!("CCF Order Idx " * string(order_idx) * " Obs " * string(obs_idx))
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
  		y = scale_order_ccf_depths ? 1. .- (1.0 .-order_ccfs_norm[:,ord_idx,obs_idx])./ccf_order_depths[ord_idx,obs_idx] :
		order_ccfs_norm[:,ord_idx,obs_idx]
		#plot!(plt2,data["v_grid"],order_ccfs_norm[:,ord_idx,obs_idx], ms=1.2,color=pal2[i],legend=:none)
		plot!(plt2,data["v_grid"],y, ms=1.2,color=pal2[i],legend=:none)
  end
  xlims!(plt2,vmin,vmax)
  xlabel!(plt2,"v (m/s)")
  ylabel!(plt2,scale_order_ccf_depths ? "Scaled CCFs" : "Normalized CCFs")
  title!(plt2,"CCF Order Idxs " * string(order_idx_lo) * "-" * string(order_idx_hi) * " Obs " * string(obs_idx))
  #ylims!(plt2,extrema(f_filtered[λmin.<=lambda.<=λmax]))
  plt2
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

# ╔═╡ 964d83d4-7c56-43b2-a989-f7d22089ee06
data["order_ccfs"][:,39,2]

# ╔═╡ c7a9d874-d823-45ff-8214-da400d01f383
collect(keys(data))

# ╔═╡ 0e4e2de1-d706-4019-bece-e24367b38d33
heatmap(data["v_grid"],obs_idx_lo:obs_idx_hi,order_ccfs_norm[:,order_idx,:].-order_ccfs_mean[:,order_idx]

# ╔═╡ ceae0ff5-f650-4dbf-81e5-26864eb1adcc
size(order_ccfs_mean)

# ╔═╡ Cell order:
# ╠═42b4e6fe-a897-11eb-2676-f7fd96f35a22
# ╠═5313b674-a897-11eb-3820-218f26a14d4d
# ╟─37a6c623-e901-4f8b-8989-835958471549
# ╠═d1fe0408-cccd-45c9-a5ef-c13a42386e44
# ╟─28aced5a-3100-4523-9f1a-f7b14016791b
# ╟─abbc4880-fedd-4bdd-a1db-b94e8a4d5654
# ╟─6c469209-c7e0-41f2-97e3-db587b894e56
# ╟─ab0fc462-a897-11eb-1727-e71788522c02
# ╟─c403d4ff-3035-4bfc-9b9e-3f33ed0a6b7a
# ╟─d84af1b1-1023-4493-8b0f-90c9030be4bd
# ╟─7dc1bf77-174b-4614-8aea-15cad25c0335
# ╟─df7571c9-8b3f-4503-8bbe-42a77f0f9c0e
# ╟─468fe2c8-7c4c-45f2-9577-a2cecba64d2d
# ╠═964d83d4-7c56-43b2-a989-f7d22089ee06
# ╠═c7a9d874-d823-45ff-8214-da400d01f383
# ╠═0e4e2de1-d706-4019-bece-e24367b38d33
# ╠═ceae0ff5-f650-4dbf-81e5-26864eb1adcc
