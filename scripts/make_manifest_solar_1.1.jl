verbose = true
 if verbose && !isdefined(Main,:RvSpectML)  println("# Loading RvSpecML")    end
 using RvSpectMLBase
 using EchelleInstruments, EchelleInstruments.NEID
 using RvSpectML
 if verbose println("# Loading NeidSolarScripts")    end
 using NeidSolarScripts
 using SunAsAStar
 using SunAsAStar.SolarRotation
 using SunAsAStar.DifferentialExtinction
# using NeidSolarScripts.SolarRotation
# using NeidSolarScripts.DifferentialExtinction
 if verbose   println("# Loading other packages")    end
 using CSV, DataFrames, Query, Dates
 #using StatsBase, Statistics, Dates

fits_target_str = "Sun"
 if !@isdefined(target_subdir)
      target_subdir = ""
 end
 if !@isdefined(create_missing_continuum_files)
     create_missing_continuum_files = false
 end
 paths_to_search_for_param = [pwd(),pkgdir(NeidSolarScripts)]

 if verbose println("# Finding what data files are avaliable.")  end
 eval(read_data_paths(paths_to_search=paths_to_search_for_param))
 @assert isdefined(Main,:neid_data_path)
 @assert isdefined(Main,:output_dir)
 output_path = joinpath(neid_data_path, target_subdir,output_dir)
 manifest_filename = joinpath(output_path,"manifest.csv")
 manifest_calib_filename = joinpath(output_path,"manifest_calib.csv")

can_skip_generating_manifest = false
if isfile(manifest_filename) || islink(manifest_filename)
    df_files  = CSV.read(manifest_filename, DataFrame)
    warn_link =  islink(manifest_filename) ? " (a link)" : ""
    println("# Read ", manifest_filename, warn_link, ".")
    @assert size(df_files,1) >= 1
    @assert hasproperty(df_files,:Filename)
    @assert hasproperty(df_files,:target)
    @assert hasproperty(df_files,:bjd)
    @assert hasproperty(df_files,:ssbz)
    @assert hasproperty(df_files,:exptime)
    @assert hasproperty(df_files,:alt_sun)
    @assert hasproperty(df_files,:Δv_diff_ext)
    @assert hasproperty(df_files,:Δfwhm²)
    #=
    @assert hasproperty(df_files,:order_snrs)
    if eltype(df_files[!,:order_snrs]) == String
        df_files[!,:order_snrs] = map(i->parse.(Float64,split(df_files[i,:order_snrs][2:end-1],',')),1:size(df_files,1))
    end
    @assert eltype(df_files[!,:order_snrs]) == Vector{Float64}
    =#
    println("# Required fields present ", size(df_files), ".  No need to regenerate manifest.")
    can_skip_generating_manifest = true
    if any(isnan.(df_files[!,:Δfwhm²]))    ||
       any(isnan.(df_files[!,:Δv_diff_ext]))   
       # any(isnan.(df_files[!,:order_snrs])) 
          can_skip_generating_manifest = false
    end
else
    println("# Will need to generate manifest at ", manifest_filename, ".")
end

#can_skip_generating_manifest = false
if can_skip_generating_manifest && !create_missing_continuum_files
    exit()
end

df_files = NEID.make_manifest(joinpath(neid_data_path, target_subdir))
df_pyrohelio_files = Pyroheliometer.make_pyrohelio_file_dataframe(pyrohelio_data_path)

calibration_target_substrings = ["Etalon","LFC","ThArS","LDLS","FlatBB"]
df_files_calib = df_files |>
        @filter( any(map(ss->contains(_.target,ss),calibration_target_substrings)) ) |>
        DataFrame

df_files_calib[!,:airmass] = fill(NaN,size(df_files_calib,1))
CSV.write(manifest_calib_filename, df_files_calib)


df_files_obs = df_files |>
    @filter( _.target == "Sun" ) |>
    @filter( !any(map(ss->contains(_.target,ss),calibration_target_substrings)) ) |>
    # @take(10) |> 
    DataFrame

#=
df_files_use = df_files_obs |>
      @filter( _.target == fits_target_str ) |>
      #@filter(bjd_first_good <= _.bjd < bjd_last_good) |>
      #@filter( is_good_day(_.Filename) ) |>
      #@take(max_spectra_to_use) |>
      DataFrame
=#

if fits_target_str == "Sun" || fits_target_str == "Solar"
    df_pyrohelio_obs = DataFrame(map(f->Pyroheliometer.get_pyrohelio_summary(df_pyrohelio_files, f), df_files_obs.Filename  ))
    df_files_obs.filename = basename.(df_files_obs.Filename)
    df_files_use = leftjoin(df_files_obs, df_pyrohelio_obs, on=:filename, makeunique=true)

    df_sol = DataFrame(get_solar_info.(df_files_use.bjd,obs=:WIYN))
    df_files_use[!,:alt_sun] = df_sol[!,:alt]
    df_files_use[!,:airmass] = df_sol[!,:airmass]
    df_files_use[!,:hour_angle] = df_sol[!,:hour_angle]
    df_files_use[!,:sol_dist] = df_sol[!,:sol_dist_au]
    #df_files_use[!,:alt_sun] = calc_solar_alt.(df_files_use.bjd,obs=:WIYN)
    #df_files_use[!,:airmass] = DifferentialExtinction.f_airmass.(DifferentialExtinction.f_z.(deg2rad.(df_files_use.alt_sun)))
    println("# Computing differential extinction")
    try
        df_files_use[!,:Δv_diff_ext] = calc_Δv_diff_extinction.(df_files_use.bjd, obs=:WIYN)
    catch
       try
          println("# Encountered error with call to JplHorizons, sleeping and will try again.")
          sleep(60)
          df_files_use[!,:Δv_diff_ext] = calc_Δv_diff_extinction.(df_files_use.bjd, obs=:WIYN)
       catch
          println("# Encountered error with call to JplHorizons, using NaNs to proceed with rest of pipeline.")
          df_files_use[!,:Δv_diff_ext] = fill(NaN,length(df_files_use.bjd))
       end
    end
    println("# FWHM effect")
    #df_files_use[!,:Δfwhm²] = CalculateFWHMDifference_SolarRotation_from_obs(df_files_use.bjd,obs=:WIYN)
    df_files_use[!,:Δfwhm²]= calc_Δfwhm_solar_rotation.(df_files_use.bjd,obs=:WIYN)

else
    @warn("This script is intended for analyzing NEID solar observations.")
    df_files_use = df_files_obs |>
      @filter( _.target == fits_target_str ) |>
      #@filter(bjd_first_good <= _.bjd < bjd_last_good) |>
      #@filter( is_good_day(_.Filename) ) |>
      #@take(max_spectra_to_use) |>
      DataFrame
end



if create_missing_continuum_files
   df_files_use.continuum_filename = map(fn->joinpath(output_path,"continuum",match(r"(neidL[1,2]_\d+[T_]\d+)\.fits$", fn)[1] * "_continuum=afs.jld2"),df_files_use.Filename)
end

#=
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


#using NaNMath
compute_order_snr = false
if compute_order_snr
println("# Computing order SNRs")
df_files_use[!,:order_snrs] = fill(zeros(0),size(df_files_use,1))
   for (i,row) in enumerate(eachrow(df_files_use))
       #if row.target != "Sun" continue   end
       if any(map(ss->contains(row.target,ss),calibration_target_substrings))
           println("# Warning: target = ", row.target)
           continue
       end
       spec = NEID.read_data(row)
       all_orders = min_order(NEID2D()):max_order(NEID2D())
       pixels = min_pixel_in_order(NEID2D()):max_pixel_in_order(NEID2D())
       order_snr = map(ord->RvSpectMLBase.calc_snr(spec, pixels, ord), all_orders)
       df_files_use[i,:order_snrs] = order_snr
   end
end

println("# Writing manifest file.")
CSV.write(manifest_filename, df_files_use)
      # Can extract vector of Reals turned into a string via:
      # df2 = CSV.read("test.csv", DataFrame)
      # map(i->parse.(Float64,split(df2[i,:order_snrs][2:end-1],',')),1:size(df2,1))


#snr_matrix = hcat(df_files_use[!,:order_snrs]...)

#=
using JLD2, FileIO
if create_missing_continuum_files
   println("# Creating missing continuum files.")
   for row in eachrow(df_files_use)
     if isfile(row.continuum_filename) continue end
     println("# Need to make ", row.continuum_filename )

     output_filename = row.continuum_filename
     tmpdir = basename(output_filename)
     #= 
     if !isdir(tmpdir)   
        println("# Making ", tmpdir)
        mkdir(tmpdir)  
     end
     =#
     #continue
     spec = NEID.read_data(row.Filename)
     try 
        continuum = Continuum.calc_continuum_model(spec)
        @save output_filename continuum
     catch
        println("# ERROR computing continuum for ", output_filename)
     end
   end
end

=#

