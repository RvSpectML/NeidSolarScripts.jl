if occursin(r"RvSpectMLEcoSystem$", pwd())
    cd("NeidSolarScripts")
    using Pkg
    Pkg.activate(".")
elseif occursin(r"NeidSolarScripts$", pwd())
   using Pkg
   Pkg.activate(".")
 elseif occursin(r"examples$", pwd())
    cd("..")
    using Pkg
    Pkg.activate(".")
 end

using JLD2, FileIO
using CSV, DataFrames, Query
#using StatsBase, Statistics
using Dates, NaNMath
using RvSpectMLBase
using EchelleInstruments

using Plots

target_subdir = "good_days/DRPv0.7"   # USER: Replace with directory of your choice
  fits_target_str = "Sun"
  output_dir = "output/continuum"
  #outputs = Dict{String,Any}()
  paths_to_search_for_param = [pwd(),"examples",joinpath(pkgdir(RvSpectMLBase),"..","RvSpectML","examples"), "/gpfs/scratch/jpn23/"]
  # NOTE: make_manifest does not update its paths_to_search when default_paths_to_search is defined here, so if you change the line above, you must also include "paths_to_search=default_paths_to_search" in the make_manifest() function call below
  pipeline_plan = PipelinePlan()
  dont_make_plot!(pipeline_plan, :movie)

verbose = false
reset_all_needs!(pipeline_plan)
#if need_to(pipeline_plan,:read_spectra)
if verbose println("# Finding what data files are avaliable.")  end
if isfile("manifest.csv")
    df_files  = CSV.read("manifest.csv", DataFrame)
    @assert size(df_files,1) >= 1
    @assert hasproperty(df_files,:Filename)
    @assert hasproperty(df_files,:target)
    @assert hasproperty(df_files,:bjd)
    @assert hasproperty(df_files,:ssbz)
    @assert hasproperty(df_files,:exptime)
else
    eval(read_data_paths(paths_to_search=paths_to_search_for_param))
    @assert isdefined(Main,:neid_data_path)
    df_files = make_manifest(neid_data_path, target_subdir, NEID )
    CSV.write("manifest.csv", df_files)
end

idx_day_to_use = 1
if verbose println("# Reading in customized parameters from param.jl.")  end
   if !@isdefined(idx_day_to_use)
       idx_day_to_use = 1
   end
   eval(code_to_include_param_jl(paths_to_search=paths_to_search_for_param))
   #=
   if match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1] ==  match(r"neidL1_(\d+)[T_](\d+)\.fits$", last(df_files_use.Filename))[1]
      match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1]date_str = match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1]
    else
      date_str = string(match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1]) * "-" * string(match(r"neidL1_(\d+)[T_](\d+)\.fits$", last(df_files_use.Filename))[1])
   end
   outputs["df_files_use"] = df_files_use

   outputs_filename = joinpath(output_dir,"solar_" * date_str * "_new.jld2")
   if isfile(outputs_filename) && false
     times_already_processed = load(outputs_filename, "times")
     files_in_day_to_process = size(df_files_solar_by_day.data[idx_day_to_use],1)
      if files_in_day_to_process == length(times_already_processed)
         println("# Already processed all ", length(times_already_processed), " files for ", date_str)
         exit()
      end
   end
   =#

continua = Vector{Array{Float32,2}}()
 for row in eachrow(df_files_use)
    m = match(r"(neidL1_\d+[T_]\d+)\.fits$", row.Filename)
    continuum_filename = joinpath(output_dir, m.captures[1] * ".jld2")
    jldopen(continuum_filename,"r") do file
        push!(continua,file["continuum"])
    end
 end

using Plots
using StatsBase

ord = 90
 pix = get_pixel_range(NEID2D(),ord)
 plt = plot()
 for obs in vcat(100:105,150:155,200:205)
    plot!(continua[obs][pix,ord], label=string(obs) )
 end
 display(plt)

 nobs = size(df_files_use,1)
 mean_continuum = mapreduce(obs->continua[obs][get_pixel_range(NEID2D(),ord),ord],.+,1:nobs) ./ nobs

 plt = plot()#legend=:none)
 mean_mean_continuum = Float64[]
 for obs in 1:nobs
    #plot!(continua[obs][pix,ord]./mean_continuum, label=string(obs) )
    push!(mean_mean_continuum,NaNMath.mean((continua[obs][pix,ord]./mean_continuum)[1000:end-1000] ))
 end
 scatter!(plt, 1:nobs,mean_mean_continuum)
 #display(plt)

 plt = plot(legend=:none)
  for obs in 1:nobs# vcat(100:105,150:155,200:205)
     pix = get_pixel_range(NEID2D(),ord)
     plot!(continua[obs][pix,ord]./(mean_continuum.*mean_mean_continuum[obs]), label=string(obs) )
  end
  display(plt)

ylims!(plt,0.99,1.01)
ylims!(plt,0.9,1.1)
