if occursin(r"RvSpectMLEcoSystem$", pwd())
    cd("NeidSolarScripts")
    using Pkg
    Pkg.activate(".")
 end

verbose = true
 make_plots = true
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectMLBase
 using EchelleInstruments, EchelleInstruments.NEID
 using EchelleCCFs
 #=
 using RvSpectML
 using NeidSolarScripts
 using NeidSolarScripts.SolarRotation
 if verbose   println("# Loading other packages")    end
 =#
 using CSV, DataFrames, Query, StatsBase, Statistics, Dates
 using JLD2, FileIO
 using NaNMath

target_subdir = "good_days/DRPv0.7"   # USER: Replace with directory of your choice
  fits_target_str = "Sun"
  output_dir = "output/continuum"
  #outputs = Dict{String,Any}()
  paths_to_search_for_param = [pwd(),"examples",joinpath(pkgdir(RvSpectMLBase),"..","RvSpectML","examples"), "/gpfs/scratch/jpn23/"]
  # NOTE: make_manifest does not update its paths_to_search when default_paths_to_search is defined here, so if you change the line above, you must also include "paths_to_search=default_paths_to_search" in the make_manifest() function call below
  pipeline_plan = PipelinePlan()
  dont_make_plot!(pipeline_plan, :movie)

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

if verbose println("# Reading in customized parameters from param.jl.")  end
   if !@isdefined(idx_day_to_use)
       idx_day_to_use = 1
   end
   eval(code_to_include_param_jl(paths_to_search=paths_to_search_for_param))
   if match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1] ==  match(r"neidL1_(\d+)[T_](\d+)\.fits$", last(df_files_use.Filename))[1]
      date_str = match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1]
    else
      date_str = string(match(r"neidL1_(\d+)[T_](\d+)\.fits$", first(df_files_use.Filename))[1]) * "-" * string(match(r"neidL1_(\d+)[T_](\d+)\.fits$", last(df_files_use.Filename))[1])
   end
   #=
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

using Distributed
addprocs(4)

@everywhere using RCall
@everywhere afs_src = joinpath(pwd(),"src","AFS.R")
@everywhere R"source($afs_src)"
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using RvSpectMLBase
@everywhere using EchelleInstruments

@everywhere function calc_continuum_model(spectrum::AbstractSpectra2D; order_idx::Integer )
    possible_pix = get_pixel_range(get_inst(spectrum),order_idx)
    bad_pix = bad_col_ranges(get_inst(spectrum),order_idx)
    pix_rng = EchelleInstruments.calc_complement_index_ranges(possible_pix,bad_pix)
    pix = mapreduce(p->collect(p),vcat,pix_rng)
    afs_inputs = zeros(Float64,length(pix),2)
    afs_inputs[:,1] .= spectrum.λ[pix,order_idx]
    afs_inputs[:,2] .= spectrum.flux[pix,order_idx]
    @assert !any(isnan.(afs_inputs))
    #=
    wv = mapreduce(p->spec.λ[p,order_idx],vcat,pix_rng)
    @assert !any(isnan.(wv))
    inten = mapreduce(p->convert(Vector{Float64},spec.flux[p,order_idx]),vcat,pix_rng)
    @assert !any(isnan.(inten))
    afs_inputs = hcat(wv,inten)
    =#
    #df = DataFrame("wv"=>wv,"intes"=>inten)
    afs_output_R = R"AFS($afs_inputs,0.95,0.25)"
    afs_output = rcopy(afs_output_R) 
    continuum = zeros(eltype(spectrum.flux),size(spectrum.flux[:,order_idx]))
    continuum = fill(NaN,size(spectrum.flux[:,order_idx]))
    continuum[pix] .= afs_output
    return continuum
end



@everywhere using EchelleInstruments.NEID
@everywhere function calc_continuum_model(spectrum::AbstractSpectra2D )
    vec_of_orders = pmap(ord->calc_continuum_model(spectrum,order_idx=ord), min_order(get_inst(spectrum)):max_order(get_inst(spectrum)) )
    output = fill(NaN, size(spectrum.flux))
    for (i,ord) in enumerate(min_order(get_inst(spectrum)):max_order(get_inst(spectrum)))
        output[:,ord] .= vec_of_orders[i]
    end
    return output
end

num_days_to_process = size(df_files_solar_by_day,1)
for idx_day_to_use in 1:num_days_to_process
  df_files_use = df_files_solar_by_day[idx_day_to_use,:data] |> @orderby(_.bjd) |> DataFrame

  println("# *** Working on day ", idx_day_to_use, " with ", size(df_files_use,1), "files. ***" )
  for row in eachrow(df_files_use)
    spec = NEID.read_data(row)
    m = match(r"neidL1_(\d+)[T_](\d+)\.fits$",row.Filename)
    output_filename = "neidL1_" * m.captures[1] * "T" * m.captures[2] * ".jld2"
    output_filename = joinpath(output_dir,output_filename)
    println("# Working on ", output_filename)
    continuum = calc_continuum_model(spec)
    @save output_filename continuum
  end

end
