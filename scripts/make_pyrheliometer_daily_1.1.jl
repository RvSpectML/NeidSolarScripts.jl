using NeidSolarScripts
#using ArgParse

#=
function parse_commandline_make_pyrheliometer_daily()
     s = ArgParseSettings( description = "Make pyrheliometer.csv for specified day.")
     @add_arg_table! s begin
         "manifest_or_date"
            help = "CSV Filename with fields l?filename, object and exptime to get pyrheliometer data for OR Date in YYYYMMDD format to extract pyrheliometer data."
            arg_type = String
            #default = "meta_test.csv" # joinpath(pwd(),"pyrohelio.csv")
            #default = "20210505"
            #default = "20210605"
         "--output"
            help = "Path for outputs: pyrheliometer.csv"
            arg_type = String
            default = joinpath(pwd(),"pyrohelio.csv")
         "--pyrheliometer_dir"
            help = "Path for pyroheliometer inputs: *.tel"
            arg_type = String
            #default = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/pyrohelio/202110_fromChad/"
            default = "/mnt/data_simons/NEID/pyrohelio/"
         "--user"
            help = "NExScI Archive username"
            arg_type = String
         "--password"
            help = "NExScI Archive password"
            arg_type = String
         "--nexsci_login_filename"
            help = "TOML file with NExScI username and password"
            arg_type = String
            default = "nexsci_id.toml"
         "--work_dir"
            help = "Working directory for NExScI query and L0 downloads."
            arg_type = String
            default = "L0"
         "--cookie"
            help = "Filename for NExScI cookie."
            arg_type = String
            default = "neidadmincookie.txt"
         "--query_filename"
            help = "Filename for L0 query results."
            arg_type = String
            default = "meta_l0.csv"
         "--root"
            help = "Path to root of data directories"
            arg_type = String
            default = ""
         "--try_tel_first"
            help = "Try pyrheliometer telemetry file before L0s"
            action = :store_true
         "--verbose"
            help = "Verbose outputs"
            action = :store_true
      end

  return parse_args(s)
end
=#
if true # verbose 
    println("# About to parse command line...") 
    flush(stdout) 
end
args = parse_commandline_make_pyrheliometer_daily()
verbose   = args["verbose"]
if verbose 
    println("# Parsed.")  
    flush(stdout) 
end
root_dir = args["root"]
output_fn = joinpath(root_dir,args["output"]) # if haskey(args,"output")     output_fn          = args["output"]      end
mkpath(dirname(output_fn))
tmp_path  = !isnothing(args["work_dir"]) ? joinpath(root_dir,args["work_dir"]) : Filesystem.mktempdir()
pyrohelio_dir = !isnothing(args["pyrheliometer_dir"]) ? joinpath(root_dir,args["pyrheliometer_dir"]) : nothing
nexsci_login_fn = !isnothing(args["nexsci_login_filename"]) ? joinpath(root_dir,args["nexsci_login_filename"]) : nothing
#verbose   = args["verbose"]

using Dates
using CSV, DataFrames, Query

if !occursin(r"\.csv$", args["manifest_or_date"])
   manifest_fn = nothing
   m = match(r"(\d{4})(\d{2})(\d{2})", args["manifest_or_date"])
   if isnothing(m)
      @error "Invalid arguement.  Must be a .csv file or date."
   end
   #println("m = ", m)
   @assert length(m.captures) == 3
   year  = parse(Int64,m[1])
   month = parse(Int64,m[2])
   day   = parse(Int64,m[3])
   @assert 2020 <= year <= 2040
   @assert 1 <= month <= 12
   @assert 1 <= day <= 31
   target_datestr = args["manifest_or_date"]
   target_date = Date(m[1]*m[2]*m[3],DateFormat("yyyymmdd"))
   df_in = DataFrame()
else
    if verbose 
        println("# About to read existing manifest...") 
        flush(stdout)
    end
   manifest_fn = joinpath(root_dir,args["manifest_or_date"])
   @assert filesize(manifest_fn) > 0
   df_in = CSV.read(manifest_fn,DataFrame)
    if verbose 
        println("# Done reading existing manifest.") 
        flush(stdout)
    end
   if hasproperty(df_in, :l0filename)   df_in = df_in |> @rename(:l0filename => :Filename) |> DataFrame  end
   if hasproperty(df_in, :l1filename)   df_in = df_in |> @rename(:l1filename => :Filename) |> DataFrame  end
   if hasproperty(df_in, :l2filename)   df_in = df_in |> @rename(:l2filename => :Filename) |> DataFrame  end
   if hasproperty(df_in, :filename)     df_in = df_in |> @rename(:filename => :Filename) |> DataFrame    end
   m = match(r"(\d{4})\/(\d{2})\/(\d{2})\/",manifest_fn)
   if isnothing(m) || (length(m.captures) != 3)
      target_datestr = nothing
      target_date = nothing
      println("# Empty ", manifest_fn, " and can't figure out date to query.  Giving up.")
      error(0)
   else
      @assert length(m.captures) == 3
      year  = parse(Int64,m[1])
      month = parse(Int64,m[2])
      day   = parse(Int64,m[3])
      @assert 2020 <= year <= 2040
      @assert 1 <= month <= 12
      @assert 1 <= day <= 31
      target_datestr = args["manifest_or_date"]
      target_date = Date(m[1]*m[2]*m[3],DateFormat("yyyymmdd"))
   end
end
if verbose 
    println("# size(df_in,1) = ", size(df_in,1)) 
    flush(stdout)  
end

using TOML
using NeidArchive

if (size(df_in,1) == 0) && (!isnothing(target_date))
   if isnothing(args["user"]) || isnothing(args["password"])
      nexsci_dict = TOML.parsefile(nexsci_login_fn)
      user_nexsci = nexsci_dict["user"]
      password_nexsci = nexsci_dict["password"]
      true
   else
      user_nexsci = args["user"]
      password_nexsci = args["password"]
      true
   end
   cookie_fn = args["cookie"]
   @assert !isnothing(user_nexsci) && !isnothing(password_nexsci) && !isnothing(cookie_fn)
   if !isdir(tmp_path)  mkpath(tmp_path) end
   query_fn  =  !isnothing(args["query_filename"]) ? joinpath(tmp_path,args["query_filename"]) : nothing

    if verbose 
        println("# About to login to Neid archive...") 
        flush(stdout)
    end
   NeidArchive.login(userid=user_nexsci, password=password_nexsci, cookiepath=cookie_fn)
    if verbose 
        println("# Done querying Neid archive.") 
        flush(stdout)
    end
   param = Dict{String,String}()
   param["datalevel"] = "solarl0"
   param["piname"] = "Mahadevan"
   param["object"] = "Sun"
   param["datetime"] = NeidArchive.datetime_one_day_solar(target_date)

    if verbose 
        println("# About to query Neid archive...") 
        flush(stdout)
    end
   NeidArchive.query(param, cookiepath=cookie_fn, outpath=query_fn)
    if verbose 
        println("# Done querying Neid archive.") 
        flush(stdout)
    end
   num_lines = countlines(query_fn) - 1
   println("# Query resulted in file with ", num_lines, " entries.")
   manifest_fn = query_fn
   @assert filesize(manifest_fn) > 0
   df_in = CSV.read(manifest_fn,DataFrame)
   df_in = df_in |> @rename(:l0filename => :Filename) |> DataFrame
end

if args["try_tel_first"] && (!isnothing(pyrohelio_dir) && isdir(pyrohelio_dir))
    if verbose 
        println("# About to make_pyrohelio_file_dataframe (1)...") 
        flush(stdout)  
    end
   pyrohelio_files = Pyroheliometer.make_pyrohelio_file_dataframe(pyrohelio_dir)
    if verbose
        println("# Done.")
        flush(stdout) 
    end
   @assert hasproperty(df_in,"Filename")
   @assert hasproperty(df_in,"exptime")
    if verbose 
        println("# About to get_pyrohelio_summary (1)...") 
        flush(stdout)  
    end
    df_out = DataFrame(map(r->Pyroheliometer.get_pyrohelio_summary(pyrohelio_files, string(r.Filename), r.exptime), eachrow(df_in)))
    if verbose
        println("# Done.")
        flush(stdout) 
    end
else
   pyrohelio_files = DataFrame()
   df_out = DataFrame()
end

if (size(df_out,1) == 0) || any(ismissing.(df_out.mean_pyroflux))
   if !isdir(tmp_path)  mkpath(tmp_path) end
   #println("# start_row = ", 1, " end_row = ", num_lines-1)
   #NeidArchive.download(query_fn, param["datalevel"], outdir=tmp_path, cookiepath=cookie_fn, start_row=1, end_row=10)
   #NeidArchive.download(query_fn, "l0", outdir=tmp_path, cookiepath=cookie_fn, start_row=1, end_row=10)
   #NeidArchive.download(query_fn, "l0", outdir=tmp_path, cookiepath=cookie_fn)
    if verbose 
        println("# About to get_pyrohelio_summary (2) into ", tmp_path, "...") 
        flush(stdout)  
    end
    df_out = DataFrame(map(fn->Pyroheliometer.get_pyrohelio_summary(joinpath(tmp_path,replace(string(fn),"neidL2_"=>"neidL0_"))),df_in.Filename))
    if verbose
        println("# Done.")
        flush(stdout) 
    end
end

if (size(df_out,1) == 0) && !args["try_tel_first"] && (!isnothing(pyrohelio_dir) && isdir(pyrohelio_dir))
    if verbose 
        println("# About to make_pyrohelio_file_dataframe (2)...") 
        flush(stdout)  
    end
   pyrohelio_files = Pyroheliometer.make_pyrohelio_file_dataframe(pyrohelio_dir)
    if verbose
        println("# Done.")
        flush(stdout) 
    end
   @assert hasproperty(df_in,"Filename")
   @assert hasproperty(df_in,"exptime")
    if verbose 
        println("# About to get_pyrohelio_summary (3)...") 
        flush(stdout)  
    end
   df_out = DataFrame(map(r->Pyroheliometer.get_pyrohelio_summary(pyrohelio_files, string(r.Filename), r.exptime), eachrow(df_in)))
    if verbose
        println("# Done.")
        flush(stdout) 
    end
end

if verbose 
    println("# About to write output...") 
    flush(stdout)
end
CSV.write(output_fn, df_out)
if verbose 
    println("# Done.")
    flush(stdout)
end
df_to_rm = df_in |> @filter( occursin(r"neidL0_\d{8}T\d+\.fits",_.Filename) ) |> @select(:Filename) |> @filter(isfile(_.Filename)) |> DataFrame
#map(fn->println("rm $fn"),df_to_rm.Filename)
exit(0)

#=

if (size(pyrohelio_files,1)>0) && (size(df_in,1)>0) # !isnothing(manifest_fn) && (filesize(manifest_fn) > 0)
      if hasproperty(df_in,"object")
         df_in = df_in |> @filter( _.object == "Sun" ) |> DataFrame
      elseif hasproperty(df_in,"target")
         df_in = df_in |> @filter( _.target == "Sun" ) |> DataFrame
      end
      filenamelist = hasproperty(df_in,"l2filename") ? df_in[!,"l2filename"] :
                     hasproperty(df_in,"l1filename") ? df_in[!,"l1filename"] :
                     hasproperty(df_in,"l0filename") ? df_in[!,"l0filename"] :
                     hasproperty(df_in,"filename") ? df_in[!,"filename"] :
                     hasproperty(df_in,"Filename") ? df_in[!,"Filename"] : nothing
      @assert !isnothing(filenamelist)
end


df_in = CSV.read(manifest_fn,DataFrame)

filenamelist = hasproperty(df_in,"l2filename") ? df_in[!,"l2filename"] :
               hasproperty(df_in,"l1filename") ? df_in[!,"l1filename"] :
               hasproperty(df_in,"l0filename") ? df_in[!,"l0filename"] :
               hasproperty(df_in,"filename") ? df_in[!,"filename"] :
               hasproperty(df_in,"Filename") ? df_in[!,"Filename"] : nothing
df_out = DataFrame(map(fn->Pyroheliometer.get_pyrohelio_summary(string(fn)),filenamelist))


using FITSIO
using Statistics
(;filename=basename(fn), time_start, exptime, mean_Δt, mean_pyroflux=mean_flux, rms_pyroflux=rms_flux, min_pyroflux=min_flux, max_pyroflux=max_flux)

df_in.mean_flux = Vector{Union{Float64,Missing}}(missing,size(df_in,1))
df_in.rms_flux = Vector{Union{Float64,Missing}}(missing,size(df_in,1))
df_in.min_flux = Vector{Union{Float64,Missing}}(missing,size(df_in,1))
df_in.max_flux = Vector{Union{Float64,Missing}}(missing,size(df_in,1))
for i in 1:size(df_in,1)
   file = hasproperty(df_in,"l2filename") ? df_in[i,"l2filename"] :
                  hasproperty(df_in,"l1filename") ? df_in[i,"l1filename"] :
                  hasproperty(df_in,"l0filename") ? df_in[i,"l0filename"] :
                  hasproperty(df_in,"filename") ? df_in[i,"filename"] :
                  hasproperty(df_in,"Filename") ? df_in[i,"Filename"] : nothing
   try
      m = match(r"neidL?_(\d+T\d+).fits",fninfull)
   	@assert m != nothing
   	time_start = DateTime(m[1],DateFormat("yyyymmddTHHMMSS"))

      f = FITS(fninfull)
      hdr = read_header(f[1])
   	exptime = hdr["EXPTIME"]
      t = read(f[21],"Time")
      v = read(f[21],"Voltage")
      fd = read(f[21],"FluxDensity")
      df_tmp = DataFrame(:Time=>t, :Voltage=>v, :FluxDensity=>fd)


      mean_Δt = get_pyrohelio_mean_Δt(df_tmp)
      #CSV.write(fnout, df)
      df_in.mean_flux[i] = mean(df_tmp.FluxDensity)
      df_in.rms_flux[i] = sqrt(var(df_tmp.FluxDensity, corrected=false))
      (df_in.min_flux[i], df_in.max_flux[i]) = extrema(df_tmp.FluxDensity)
   catch ex
      println("# Failed to extract pyrheliometer data from $fninfull.\n")
      continue
  end # try to read L0 fits file
end # for each L0 file


try
   mean_flux = mean(df.FluxDensity)
   rms_flux = sqrt(var(df.FluxDensity, corrected=false))
   outputline = fnin * "," * string(mean_flux) * "," * string(rms_flux) * "\n"
   write(exampleFileIOStream, outputline)
   if mod(i,10) == 0 flush(exampleFileIOStream)  end
 catch ex
    println("# Failed to write summary data for $fn_out\n.")
    continue
 end






target_date = !isnothing(args["date"]) ? DateTime(args["date"],DateFormat("yyyymmdd")) : nothing






Pyroheliometer.pick_pyrohelio_file(pyrohelio_files,target_date)
Pyroheliometer.get_pyrohelio_summary(pyrohelio_files, fn, exptime)


if !isnothing(args["pyrheliometer_dir"])
   pyrohelio_dir = arg["pyrheliometer_dir"]
   global df_pyrohelio_files = Pyroheliometer.make_pyrohelio_file_dataframe(pyrohelio_data_path)
   # Extract from pyrheliometer telimetry file
end




using CSV, DataFrames
using NeidArchive




if haskey(args,"root")       root_dir            = args["root"]   end
if haskey(args,"input")      input_dir           = args["input"]   end
if haskey(args,"output")     output_dir          = args["output"]      end
if haskey(args,"subdir")     subdir              = args["subdir"]      end
if haskey(args,"pyrohelio")  pyrohelio_dir       = args["pyrohelio"]   end
create_missing_continuum_files = haskey(args,"create-continuum") ? args["create-continuum"] : false
verbose = haskey(args,"verbose") ? args["verbose"] : false

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
#=
 if !@isdefined(subdir)
      subdir = ""
 end
 if !@isdefined(create_missing_continuum_files)
     create_missing_continuum_files = false
 end
=#
 paths_to_search_for_param = [pwd(),pkgdir(NeidSolarScripts)]

 if verbose println("# Finding what data files are avaliable.")  end
 eval(read_data_paths(paths_to_search=paths_to_search_for_param))
 @assert isdefined(Main,:root_dir)
# @assert isdefined(Main,:neid_data_path)
 @assert isdefined(Main,:input_dir)
 @assert isdefined(Main,:output_dir)
 @assert isdefined(Main,:subdir)
 input_path = joinpath(root_dir,input_dir, subdir)
 output_path = joinpath(root_dir,output_dir, subdir)
 manifest_filename = joinpath(output_path,"manifest.csv")
 manifest_calib_filename = joinpath(output_path,"manifest_calib.csv")
 pyrohelio_data_path = joinpath(root_dir,pyrohelio_dir)

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

try
   global df_pyrohelio_files = Pyroheliometer.make_pyrohelio_file_dataframe(pyrohelio_data_path)
catch ex
   println("# Error making pyroheliometer dataframe for " * pyrohelio_data_path * ".")
   exit(0)
end

try
   global df_files = NEID.make_manifest(input_path) #joinpath(root_path,input_dir, subdir))

catch ex
   println("# Error making manifest for ", subdir," and/or pyroheliometer dataframe.")
   touch(manifest_filename)
   touch(manifest_calib_filename)
   exit(0)
end

if size(df_files,1) == 0
   println("# No files in manifest.")
   touch(manifest_filename)
   touch(manifest_calib_filename)
   exit(0)
end

calibration_target_substrings = ["Etalon","LFC","ThArS","LDLS","FlatBB"]
df_files_calib = df_files |>
        @filter( any(map(ss->contains(_.target,ss),calibration_target_substrings)) ) |>
        DataFrame

df_files_calib[!,:airmass] = fill(NaN,size(df_files_calib,1))
CSV.write(manifest_calib_filename, df_files_calib)


try
    global df_files_obs = df_files |>
    @filter( _.target == "Sun" ) |>
    @filter( !isnothing(_.drpversion) && !ismissing(_.drpversion) && !(_.drpversion == "") ) |>
    #@filter( Base.thisminor(VersionNumber(_.drpversion)) == Base.thisminor(VersionNumber(1,1,0)) ) |>
    @filter( !any(map(ss->contains(_.target,ss),calibration_target_substrings)) ) |>
    # @take(10) |>
    DataFrame
catch ex
   println("# No files passed target and DRP version criteria.")
   touch(manifest_filename)
   exit(0)
end


#=
df_files_use = df_files_obs |>
      @filter( _.target == fits_target_str ) |>
      #@filter(bjd_first_good <= _.bjd < bjd_last_good) |>
      #@filter( is_good_day(_.Filename) ) |>
      #@take(max_spectra_to_use) |>
      DataFrame
=#

if size(df_files_obs,1) == 0
   touch(manifest_filename)
   println("# No files avaliable in df_files_obs.")
   exit(0)
end

if fits_target_str == "Sun" || fits_target_str == "Solar"
    df_pyrohelio_obs = DataFrame(map(x->Pyroheliometer.get_pyrohelio_summary(df_pyrohelio_files, x[1], x[2]), zip(df_files_obs.Filename, df_files_obs.exptime)  ))
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

if size(df_files_use,1) == 0
   println("# No files avaliable.")
   touch(manifest_filename)
   exit(0)
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
        mkpath(tmpdir)
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
=#
