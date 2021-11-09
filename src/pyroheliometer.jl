module Pyroheliometer
using Dates
using DataFrames, CSV
using FITSIO

function make_pyrohelio_file_dataframe(pyro_path::String)
	pyro_files = readdir(pyro_path)
	filter!(f->match(r"neid_ljpyrohelio_chv0_(\d+).tel",f)!=nothing , pyro_files )
	df_pyrohelio_files = DataFrame(:filename=>joinpath.(pyro_path,pyro_files), :start_date=>map(f->Date(match(r"neid_ljpyrohelio_chv0_(\d+).tel$",f)[1],DateFormat("yyyymmdd")),pyro_files))
end

begin
	function pick_pyrohelio_file(pyrohelio_files::DataFrame, fn::String)
		m = match(r"neidL2_(\d+T\d+).fits",fn)
		@assert m != nothing
		t_start = DateTime(m[1],DateFormat("yyyymmddTHHMMSS"))
		pick_pyrohelio_file(pyrohelio_files,t_start)
	end

	function pick_pyrohelio_file(pyrohelio_files::DataFrame, t::DateTime)
		if !issorted(df_pyrohelio_files, :start_date)
			sort!(df_pyrohelio_files, :start_date)
		end
		idx = searchsortedlast(pyrohelio_files.start_date, t)
		if !(1 <= idx <= size(pyrohelio_files,1))
			throw("pick_pyrohelio_file: idx = " * string(idx))
		end
		return pyrohelio_files.filename[idx]
	end
end

# Cache so don't reread CSV files every time
pyrohelio_cache = Dict{String,DataFrame}()

function clear_pyrohelio_cache!()
	global pyrohelio_cache = Dict{String,DataFrame}()
end

begin
	function read_pyrohelio_file(fn::String; force::Bool = false)
		if !force && haskey(pyrohelio_cache,fn)
			df_pyrohelio_data = pyrohelio_cache[fn]
		else
			df_pyrohelio_data = CSV.read(fn, header=[:time_str,:voltage,:flux],DataFrame)

			datefmt = DateFormat("yyyy-mm-ddTHH:MM:SS.sss")
			df_pyrohelio_data.time = map(tstr->DateTime(view(tstr,1:23),datefmt),collect(df_pyrohelio_data.time_str))
			#delete!(pyrohelio_data_from_prev_file,[:time_str])
			pyrohelio_cache[fn] = df_pyrohelio_data
		end
		df_pyrohelio_data
	end

	function read_pyrohelio_file(pyrohelio_files::DataFrame, t::DateTime)
		fn = pick_pyrohelio_file(pyrohelio_files, t)
		#read_pyrohelio_file(joinpath(pyro_path,fn))
		read_pyrohelio_file(fn)
	end
	function read_pyrohelio_file(pyrohelio_files::DataFrame, fn::String)
		m = match(r"neidL2_(\d+T\d+).fits",fn)
		@assert m != nothing
		t_start = DateTime(m[1],DateFormat("yyyymmddTHHMMSS"))
		fn = pick_pyrohelio_file(pyrohelio_files, t_start)
		#read_pyrohelio_file(joinpath(pyro_path,fn))
		read_pyrohelio_file(fn)
	end
end

begin
	function get_pyrohelio_data(pyrohelio_files::DataFrame, t_start::DateTime, exptime::Real)
		t_stop = t_start + Dates.Second(exptime)

		df = read_pyrohelio_file(pyrohelio_files,t_start)
		idx_start = searchsortedfirst(df.time,t_start)
		if !(1 <= idx_start <= size(df,1))
			throw("get_pyrohelio_data: idx_start = " * string(idx_start))
		end
		idx_stop = searchsortedlast(df.time,t_stop)
		if !(1 <= idx_stop <= size(df,1))
			throw("get_pyrohelio_data: idx_stop = " * string(idx_stop))
		end
		#delta_idx = searchsortedlast(view(pyrohelio_data_from_prev_file.time,idx_start:size(pyrohelio_data_from_prev_file,1)),t_stop)
		#idx_stop = idx_start+delta_idx-1
		df[idx_start:idx_stop,:]
	end

	function get_pyrohelio_data(pyrohelio_files::DataFrame, fn::String)
		m = match(r"neidL2_(\d+T\d+).fits",fn)
		@assert m != nothing
		t_start = DateTime(m[1],DateFormat("yyyymmddTHHMMSS"))
		fits = FITS(fn)
		@assert length(fits) >= 1
		hdr = read_header(fits[1])
		exptime = hdr["EXPTIME"]
		exptime_sec = round(Int64,exptime)
		close(fits)
		get_pyrohelio_data(pyrohelio_files,t_start, exptime_sec)
	end
end

function get_pyrohelio_mean_Δt(df::DataFrame; time_start::DateTime = first(df.time))
	dt = map(dt->dt.value/1000,df.time .-time_start)
	sum(dt.*df.flux)/sum(df.flux)
end

function get_pyrohelio_summary(pyrohelio_files::DataFrame, fn::String)
	m = match(r"neidL2_(\d+T\d+).fits",fn)
	@assert m != nothing
	time_start = DateTime(m[1],DateFormat("yyyymmddTHHMMSS"))
	fits = FITS(fn)
	@assert length(fits) >= 1
	hdr = read_header(fits[1])
	exptime = hdr["EXPTIME"]
	exptime_sec = round(Int64,exptime)
	close(fits)

	try
		df = get_pyrohelio_data(pyrohelio_files, time_start, exptime_sec)
		#time_start = first(df.time)  # Comment out so use FITS file start as reference time
		mean_Δt = get_pyrohelio_mean_Δt(df, time_start=time_start)
		mean_flux = mean(df.flux)
		rms_flux = sqrt(var(df.flux,mean=mean_flux,corrected=false))
		(;filename=basename(fn), time_start, exptime, mean_Δt, mean_pyroflux=mean_flux, rms_pyroflux=rms_flux)
	catch exepction
		(;filename=basename(fn), time_start, exptime, mean_Δt=missing, mean_pyroflux=missing, rms_pyroflux=missing)
	end
end

export make_pyrohelio_file_dataframe, pick_pyrohelio_file, read_pyrohelio_file
export get_pyrohelio_data, get_pyrohelio_summary
end # module
