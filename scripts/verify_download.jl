using Markdown
using InteractiveUtils
using CSV, DataFrames, Query
using Dates, MD5
using ArgParse

  function parse_commandline_verify_downloads()
     s = ArgParseSettings( description = "Verify FITS file downloads match contents of meta.csv.")
     @add_arg_table! s begin
         "inputdir"
            help = "Directory with downloaded FITS files."
            arg_type = String
            required = true
         "--object"
            help = "Filter for object name matching (default: empty, so no checking)"
            arg_type = String
            default = ""
#=
         "--root"
             help = "Path to root of data directories"
             arg_type = String
             default = "/gpfs/group/ebf11/default/RISE_NEID_ebf/data"
         "--input"
             help = "Path to inputs FITS files (fits)"
             arg_type = String
             default = "solar_L2/v1.1"
         "--output"
             help = "Path to daily CCFs (jld2)"
             arg_type = String
             default = "output"
=#
         "--checksums"
            help = "Compute md5 checksums for input FITS files."
            action = :store_true
         "--max_num_spectra"
            help = "Maximum number of spectra to process."
            arg_type = Int
            default = 300  
         "--crawl"
            help = "Run on each date's subdirectory."
            action = :store_true
         "--quiet"
            help = "Don't write output to terminal"
            action = :store_true
         "--overrwrite"
            help = "Specify it's ok to overwrite files."
            action = :store_true
      end

     return parse_args(s)
 end
 args = parse_commandline_verify_downloads()


compute_checksums = args["checksums"]
input_path = abspath(args["inputdir"])

begin
	"""
	make_meta_redownload( input_path)
	Compare MD5 checksums of data in input_path in meta.csv from NEID archive to actual (or cached values in md5.csv)
	"""
	function make_meta_redownload end
	
	function make_meta_redownload(input_path::String ; write_md5s::Bool = true, write_meta_redownload::Bool = true, max_files::Integer = 300, checksums::Bool = true )
		meta_filename = joinpath(input_path,"meta.csv") 
		manifest_filename = joinpath(input_path,"md5.csv") 
		meta_out_filename = joinpath(input_path,"meta_redownload.csv")
		verified_filename = joinpath(input_path,"0_download_verified")
                mixed_versions_filename = joinpath(input_path,"0_mixed_versions")
		
		if isfile(verified_filename)
                        time_verified = mtime(verified_filename)
                else
                        time_verified = 0
		end
		if isfile(meta_filename)
                        time_meta = mtime(meta_filename)
                else
			@warn("Can't open $meta_filename.")
			return DataFrame()
		end
	        manifest = DataFrame(:Filename=>readdir(input_path,join=true)) |> @filter( occursin(r"\.fits$",_.Filename) ) |> DataFrame
                manifest.mtime = mtime.(manifest.Filename)
                max_fits_mtime = maximum(manifest.mtime)

		if isfile(manifest_filename) && (mtime(manifest_filename) > max_fits_mtime)
                        args["quiet"] || println("# Reusing $manifest_filename.")
                        time_manifest = mtime(manifest_filename)
                        manifest = CSV.read(manifest_filename,DataFrame)
                else
                        # In principle, could only update md5 sums for files modified since creating last md5.csv
                        time_manifest = datetime2unix(now())
                        if checksums
			   args["quiet"] || println("# Computing md5sums's for ", size(manifest,1), " files.")
                           manifest = add_md5s!(manifest, max_files=300)
                           if write_md5s
			      args["quiet"] || println("# Writing $manifest_filename.")
                              CSV.write(manifest_filename,manifest)
                           end
                        end
		end
	
                if (time_verified >= time_meta) && (time_verified >= time_manifest) && (time_verified >= max_fits_mtime)
			return DataFrame()
                end
		
                meta = CSV.read(meta_filename,DataFrame)
		#if args["target"] == "Sun"  # TODO: Eventually need to generalize for non-solar observations
	        if length(args["object"])>0 
                   meta = meta |> @filter( occursin(args["object"],_.object) ) |> DataFrame
                end
                if in("swversion" , names(meta) )
		   m = match(r"(v\d+\.\d+)\.\d+\/", input_path)
                   if m != nothing
                      drp_latest_minor = Base.thisminor(VersionNumber(m[1]))
                   else
                      drp_latest_minor = Base.thisminor(maximum(VersionNumber.(meta.swversion)))
                   end
                   if !args["quiet"] 
                      println("# Checking minor version matches: ", drp_latest_minor)
                   end
                   if ! all(Base.thisminor.(VersionNumber.(meta.swversion)) .== drp_latest_minor) 
                      msg = "meta.csv contains multiple DRP versions: " * string(unique(meta.swversion))
                      @warn(msg)
                      #touch(mixed_versions_filename)
                      open(mixed_versions_filename, "w") do io
                          write(io, msg)
                      end;
                      # TODO: Should we automatically download new metadata here?
                      # Hard since we won't know all the right query flags (perhaps save them to the directory?)
                      # For now let a script find this file and trigger a reprocess
                   end
		end
                manifest.filename = basename.(manifest.Filename)
		need_to_redownload_df = make_meta_redownload(meta,manifest, input_path, checksums=checksums,max_files=max_files)
		
		if checksums && size(need_to_redownload_df,1) == 0
			touch(verified_filename)
			rm(meta_out_filename,force=true)
		end
		if write_meta_redownload && (size(need_to_redownload_df,1)>=1)
			CSV.write(meta_out_filename,need_to_redownload_df)
		end
		need_to_redownload_df
	end
	
	function add_md5s!(manifest_df::DataFrame; max_files::Integer = 300, checksums::Bool = true)
		download_stats_df = manifest_df |> @take(max_files) |> @map({_.Filename, mtime_download=mtime(_.Filename), md5_download=string(bytes2hex(open(md5,_.Filename))) })  |> DataFrame
        end

	function make_meta_redownload(meta_df::DataFrame, manifest_df::DataFrame, path::String; max_files::Integer = 300, checksums::Bool = true, suspect_dirname = "suspect" )
                if size(manifest_df,1) == 0
                   return meta_df
                   @warn("Empty dataframe meta: " * string(size(meta_df,1)) * ", manifest: " * string(size(manifest_df,1)) )
                end
                if size(meta_df,1) == 0 || size(manifest_df,1) == 0
                   @warn("Empty dataframe meta: " * string(size(meta_df,1)) * ", manifest: " * string(size(manifest_df,1)) )
                   return DataFrame()
                end 
		if in("l2filename",names(meta_df))
		   meta_missing_files = meta_df |> @filter( !(_.l2filename ∈ basename.(manifest_df.filename)) ) |> DataFrame
		elseif in("l1filename",names(meta_df))
		   meta_missing_files = meta_df |> @filter( !(_.l1filename ∈ basename.(manifest_df.filename)) ) |> DataFrame
		else
		   meta_missing_files = meta_df |> @filter( !(_.l0filename ∈ basename.(manifest_df.filename)) ) |> DataFrame
		end
		need_to_redownload_df = copy(meta_missing_files)
		if checksums
		     if in("l2filename",names(meta_df))
                        download_success_df = manifest_df |> @join(meta_df, String(_.filename), String(_.l2filename), {_.filename, success= _.md5_download == __.l2checksum}) |> DataFrame
			meta_bad_checksum_df = DataFrame(download_success_df |> @filter(!_.success) |> @join(meta_df, String(_.filename), String(_.l2filename), {_.filename, meta=__}) |> DataFrame).meta
		     elseif in("l1filename",names(meta_df))
                        download_success_df = manifest_df |> @join(meta_df, String(_.filename), String(_.l1filename), {_.filename, success= _.md5_download == __.l1checksum}) |> DataFrame
			meta_bad_checksum_df = DataFrame(download_success_df |> @filter(!_.success) |> @join(meta_df, String(_.filename), String(_.l1filename), {_.filename, meta=__}) |> DataFrame).meta
		     else
                        download_success_df = manifest_df |> @join(meta_df, String(_.filename), String(_.l0filename), {_.filename, success= _.md5_download == __.l0checksum}) |> DataFrame
			meta_bad_checksum_df = DataFrame(download_success_df |> @filter(!_.success) |> @join(meta_df, String(_.filename), String(_.l0filename), {_.filename, meta=__}) |> DataFrame).meta
		     end
                        if size(meta_bad_checksum_df,1) >= 1 
			   append!(need_to_redownload_df,meta_bad_checksum_df)
                           suspect_files_dir = joinpath(path,suspect_dirname)
                           isdir(suspect_files_dir) || mkdir(suspect_files_dir)
                           for file in meta_bad_checksum_df.filename
                               mv(file, joinpath(suspect_files_dir,file), force=true)
                           end    
                        end 
                 end
                 begin 
                        meta_old_drp_df = DataFrame()
                        if in("swversion",names(meta_df))
			   drp_latest_minor = Base.thisminor(maximum(VersionNumber.(meta_df.swversion)))
                           meta_old_drp_df = meta_df |> @filter( Base.thisminor(VersionNumber(_.swversion)) < drp_latest_minor ) |> DataFrame
			end
                        if size(meta_old_drp_df,1) >= 1 
			   append!(need_to_redownload_df,meta_old_drp_df)
                           suspect_files_dir = joinpath(path,suspect_dirname)
                           isdir(suspect_files_dir) || mkdir(suspect_files_dir)
                           if in("swversion",names(meta_df))
                               for file in meta_old_drp_df.filename
                                   mv(file, joinpath(suspect_files_dir,file), force=true)
                               end    
			   else
                               for file in meta_old_drp_df.l0filename
                                   mv(joinpath(path,file), joinpath(suspect_files_dir,file), force=true)
                               end    
			   end
                        end
		end
		need_to_redownload_df
	end
end

need_to_redownload_df = DataFrame()

if !args["crawl"]

   if !args["quiet"] 
      println("# calling make_meta_redownload(",input_path,", ", "checksums=",compute_checksums,")" )
   end
   need_to_redownload_df = make_meta_redownload(input_path,checksums=compute_checksums)

   if !args["quiet"] 
      println("# verify_downloads.jl $input_path identified ", size(need_to_redownload_df,1), " files to download.")
   end

else
  subdirs = DataFrame(:Filename=>readdir(input_path,join=true))
  #println(subdirs.Filename)

  subdirs = subdirs |> @filter( occursin(r"\d{4}\-\d{2}\-\d{2}",_.Filename) ) |> DataFrame
  #println("# Filtered")
  #println(subdirs)

  for dir in subdirs.Filename

      if !args["quiet"] 
         println("# calling make_meta_redownload(",dir,", ", "checksums=",compute_checksums,")" )
      end
      global need_to_redownload_df = make_meta_redownload(dir,checksums=compute_checksums)
      if !args["quiet"] 
         println("# verify_downloads.jl $dir identified ", size(need_to_redownload_df,1), " files to download.")
      end
  #break
  end

end

