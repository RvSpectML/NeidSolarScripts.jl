if occursin(r"RvSpectMLEcoSystem$", pwd())   cd("NeidSolarScripts")   end
 using Pkg
 Pkg.activate(".")
 using HDF5, JLD2, FITSIO, FileIO, DataFrames, Query, RvSpectMLBase, Statistics, MultivariateStats, Plots

if occursin(r"NeidSolarScripts$", pwd())   cd("output56")   end
files = DataFrame(:filename=>readdir())
jld_files = files |> @filter(contains(_.filename,"new.jld2")) |> DataFrame

plt = plot()
 for day in eachrow(jld_files)
    #f = h5open("solar_20210104_new.jld2")
    f = h5open(day.filename)
    ccfs = read(f,"ccfs_espresso")
    M = fit(PCA,ccfs[401:end-401,60:end-60],maxoutdim=10,pratio=1.0)
    plot!(plt,M.mean./maximum(M.mean), label="Mean")
    map(i->plot!(plt,sign(M.proj[200,i]).*M.proj[:,i]./maximum(abs.(extrema(M.proj[:,i]))), label=string(i)),1:4);
    break
 end
 display(plt)
 #savefig("ccf_pca_20120104.png")

v_grid =read(h5open(jld_files.filename[1]),"v_grid")
v_grid = 100.0 * (-616:616)
jld_files.mean_ccf = Vector{Vector{Float64}}(undef,size(jld_files,1))
jld_files.smooth_ccf = Vector{Vector{Float64}}(undef,size(jld_files,1))
jld_files.ccf_resid = Vector{Array{Float64,2}}(undef,size(jld_files,1))
#plt = plot()
for (i,day) in enumerate(eachrow(jld_files))
     #f = h5open("solar_20210104_new.jld2")
     f = h5open(day.filename)
     #jld_files.mean_ccf[i] = read(f,"mean_ccf")
     #jld_files.smooth_ccf[i] = read(f,"ccf_template_smooth")
     jld_files.ccf_resid[i] = read(f,"ccf_resid_minus_rv_proj")
  end
  #display(plt)
  #savefig("ccf_pca_20120104.png")
mean_ccf_matrix = reduce(hcat,jld_files.mean_ccf)
smooth_ccf_matrix = reduce(hcat,jld_files.smooth_ccf)

norm_mean = vec((mean(mean_ccf_matrix[1:300,:],dims=1).+mean(mean_ccf_matrix[end-300:end,:],dims=1))/2)
norm_smooth = vec((mean(smooth_ccf_matrix[1:300,:],dims=1).+mean(smooth_ccf_matrix[end-300:end,:],dims=1))/2)
norm_mean = vec((mean(mean_ccf_matrix[600:632,:],dims=1))/2)
norm_smooth = vec((mean(smooth_ccf_matrix[600:632,:],dims=1))/2)

plot(v_grid,mean_ccf_matrix./norm_mean'.-smooth_ccf_matrix./norm_smooth')
plot(v_grid,(mean_ccf_matrix.-smooth_ccf_matrix)./norm_smooth')
xlims!(-12e3,12e3)

plot(v_grid,(mean_ccf_matrix-smooth_ccf_matrix)[:,1])

cols_use = 301:(size(mean_ccf_matrix,1)-301)
M = fit(PCA,mean_ccf_matrix[cols_use,:],maxoutdim=10,pratio=1.0)
principalvars(M)
scatter(log10.(principalvars(M)./tprincipalvar(M)));
 xlabel!("Number of PCs");
 ylabel!("log frac variance remaining")

plot(v_grid[cols_use,:],mean(mean_ccf_matrix[cols_use,:],dims=2),label="Mean");
 plot!(v_grid[cols_use,:],M.proj[:,1].+0.5,label="1");
 plot!(v_grid[cols_use,:],M.proj[:,2].+0.4,label="2");
 plot!(v_grid[cols_use,:],M.proj[:,3].+0.25,label="3");
 plot!(v_grid[cols_use,:],M.proj[:,4].+0.15,label="4");
 plot!(v_grid[cols_use,:],M.proj[:,5].+0.,label="5");
 xlims!(-1e4,1e4)
#savefig("daily_ccf_PCA.png")

f8 = h5open(jld_files.filename[8])
ccfs =read(f8,"ccfs_espresso")

plot(v_grid,ccfs[:,100:120].-mean(ccfs[:,100:120],dims=2))

plot(v_grid,ccfs[:,100:120].-mean(ccfs[:,100:120],dims=2))
xlims!(-1.2e3,1.2e3)
