### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 75225d8a-ad5e-11eb-131a-b76726982c01
using CSV, DataFrames, Query, Statistics, NaNMath

# ╔═╡ 91525f7c-027d-41b2-a2eb-33f9e1eddf46
using Plots

# ╔═╡ 630e63c6-d61e-46ad-b895-05db5902ad27
begin 
	data = Dict{String,Any}()
	path = "/mnt/data_simons/NEID/DRPv0.7-fixedflatfielding2/output/manifests"
	for (root, dirs, files) in walkdir(path)
	    #=
		println("Directories in $root")
		for dir in dirs
	        println(joinpath(root, dir)) # path to directories
	    end
		=#
	    println("Files in $root")
	    for file in files
			m = match(r"(\d+)/manifest\.csv",joinpath(root, file))
			if m == nothing continue end
			datestr = m.captures[1]
	        #println(joinpath(root, file)) # path to files
			println(datestr) # path to files
			df = CSV.read(joinpath(root, file),DataFrame)
			data[datestr] = df
	    end
		#break
	end
end

# ╔═╡ faf36f02-c231-4cb0-a095-d981762a2f99
begin 
	order_to_use = 60
	for (i,k) in enumerate(sort(collect(keys(data))))
		nobs = length(data[k].order_snrs)
		if nobs == 0
			pop!(data,k)
		end
	end
	ndays = length(keys(data))
	days_snrs = fill(Vector{Float64}(undef,0),ndays)
	days_airmass = fill(Vector{Float64}(undef,0),ndays)
	days_sha = fill(Vector{Float64}(undef,0),ndays)
	for (i,k) in enumerate(sort(collect(keys(data))))
		nobs = length(data[k].order_snrs)	
		days_snrs[i] = map(i->parse.(Float64,split(data[k].order_snrs[i][2:end-1],','))[order_to_use],1:length(data[k].order_snrs))
		days_airmass[i] = map(i->data[k].airmass[i],1:length(data[k].airmass))
		days_sha[i] = map(i->data[k].bjd[i].-floor(data[k].bjd[i]),1:length(data[k].bjd))
	end
end

# ╔═╡ f6545aab-0865-44aa-973a-644101b65e20
scatter(mean.(days_snrs))

# ╔═╡ 10da0778-c51c-4f98-90e6-b1d647274b98
scatter(maximum.(days_snrs))

# ╔═╡ aa3aa153-51e1-4649-af89-07246ed0c00d
scatter(std.(days_snrs))

# ╔═╡ 62e42b0f-07e9-47b7-bc2c-d52c86c65e38
threshold = maximum(maximum.(days_snrs))

# ╔═╡ 755fbea9-adcf-4c24-b837-118aaba14ff6
begin
	plt = plot(legend=:none)
	for obsid in 1:ndays
		idx_keep = days_snrs[obsid] .>= 0.80*threshold 
		scatter!(plt,days_sha[obsid][idx_keep].*24,days_snrs[obsid][idx_keep]./threshold, ms=1)
		idx_keep = days_snrs[obsid] .>= 0.92*threshold 
		scatter!(plt,days_sha[obsid][idx_keep].*24,days_snrs[obsid][idx_keep]./threshold, ms=2)
	end
	plt
end

# ╔═╡ 9db1b522-f069-4319-8563-762e33db910f
4/24

# ╔═╡ Cell order:
# ╠═75225d8a-ad5e-11eb-131a-b76726982c01
# ╠═630e63c6-d61e-46ad-b895-05db5902ad27
# ╠═faf36f02-c231-4cb0-a095-d981762a2f99
# ╠═f6545aab-0865-44aa-973a-644101b65e20
# ╠═10da0778-c51c-4f98-90e6-b1d647274b98
# ╠═aa3aa153-51e1-4649-af89-07246ed0c00d
# ╠═91525f7c-027d-41b2-a2eb-33f9e1eddf46
# ╠═755fbea9-adcf-4c24-b837-118aaba14ff6
# ╠═62e42b0f-07e9-47b7-bc2c-d52c86c65e38
# ╠═9db1b522-f069-4319-8563-762e33db910f
