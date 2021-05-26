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
	using PlutoUI
	using Plots
end

# ╔═╡ 5313b674-a897-11eb-3820-218f26a14d4d
begin
	fits_path = "/mnt/data_simons/NEID/DRPv0.7-fixedflatfielding2/"
	fits_fn = "neidL1_20210104T182500.fits"
	continuum_afs_path = joinpath(fits_path,"output","continuum")
	continuum_afs_fn = "neidL1_20210104T182500_continuum=afs.jld2"

	spec = NEID.read_data(joinpath(fits_path,fits_fn))
end;

# ╔═╡ 28aced5a-3100-4523-9f1a-f7b14016791b
md"Order index  $(@bind order_idx NumberField(1:100; default=60))
"

# ╔═╡ 8352fbc6-483c-41fa-8649-1d3d40a33ef6
md"   Stretch factor $(@bind stretch_fac NumberField(0.05:20; default=1.0))
   Merging threshold $(@bind merge_thresh NumberField(0:1; default=0.05))
   Ready $(@bind calc_continuum CheckBox(default=false))
"

# ╔═╡ f36dee02-9a27-4a58-938f-5bad6a9ba286
md"Stretch Factor = $stretch_fac   Merging Factor = $merge_thresh"

# ╔═╡ 37b20ff3-e93a-4ce1-82eb-29c23581add3
if 1<=order_idx<=100
	 lambda = spec.λ[:,order_idx]
	 f_obs = spec.flux[:,order_idx]
	 lambda_plt = lambda
	 λmin_default = minimum(lambda)
	 λmax_default = maximum(lambda)
end;

# ╔═╡ 851d39ba-a897-11eb-28e5-1191c9c9e172
if calc_continuum && 0<stretch_fac<=40 && 0 <= merge_thresh <10
	 (anchors, continuum, f_filtered_internal) = Continuum.calc_continuum(lambda,f_obs,
	      stretch_factor=stretch_fac, merging_threshold=merge_thresh, verbose=true)
	 println("# num_anchors = ", length(anchors))
	 #f_filtered = Continuum.calc_rolling_median(lambda,f_obs, width=4)
end;

# ╔═╡ f99f1957-0799-4296-a681-874b7407bb58
if calc_continuum
	f_filtered = f_filtered_internal
	continuum_plt = continuum
	anchors_plt = anchors
else
	f_filtered = Continuum.calc_rolling_median(lambda,f_obs, width=4)
	continuum_plt = nothing
	anchors_plt = nothing
end;

# ╔═╡ abbc4880-fedd-4bdd-a1db-b94e8a4d5654
@bind λmin NumberField(4000:10000; default=λmin_default)

# ╔═╡ c12d7cda-68be-4b78-b3e0-26a93cb6ac39
@bind λmax NumberField(4000:10000; default=λmax_default)

# ╔═╡ ce053df6-a899-11eb-2024-a9d147acc990
#= md"λ_min  $(@bind λmin NumberField(4000:10000; default=λmin_default))
   λ_max  $(@bind λmax NumberField(4000:10000; default=λmax_default))
"
=#

# ╔═╡ ab0fc462-a897-11eb-1727-e71788522c02
if 10000>λmax>λmin>4000
  plt = scatter(lambda,f_obs, ms=1.2, color=:grey,legend=:none)
	if f_filtered != nothing && anchors_plt != nothing
  		plot!(lambda_plt,f_filtered,color=:cyan)
  		plot!(lambda,continuum_plt,color=:magenta)
  		scatter!(lambda[anchors_plt], f_filtered[anchors_plt], color=:blue, ms=4)
	end
  #ylims!(0.95,1.15)
  #xmin = 5365
  #xmax = xmin+35
  xlims!(plt,λmin,λmax)
  ylims!(plt,extrema(f_filtered[λmin.<=lambda.<=λmax]))
  plt
end

# ╔═╡ Cell order:
# ╠═42b4e6fe-a897-11eb-2676-f7fd96f35a22
# ╠═5313b674-a897-11eb-3820-218f26a14d4d
# ╟─28aced5a-3100-4523-9f1a-f7b14016791b
# ╠═8352fbc6-483c-41fa-8649-1d3d40a33ef6
# ╟─f36dee02-9a27-4a58-938f-5bad6a9ba286
# ╟─37b20ff3-e93a-4ce1-82eb-29c23581add3
# ╟─851d39ba-a897-11eb-28e5-1191c9c9e172
# ╟─f99f1957-0799-4296-a681-874b7407bb58
# ╠═abbc4880-fedd-4bdd-a1db-b94e8a4d5654
# ╠═c12d7cda-68be-4b78-b3e0-26a93cb6ac39
# ╟─ce053df6-a899-11eb-2024-a9d147acc990
# ╠═ab0fc462-a897-11eb-1727-e71788522c02
