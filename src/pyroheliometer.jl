module Pyroheliometer
using Dates
using DataFrames, CSV
using FITSIO
using Statistics

function make_pyrohelio_file_dataframe(pyro_path::AbstractString)
	pyro_files = readdir(pyro_path)
	filter!(f->match(r"neid_ljpyrohelio_chv0_(\d+).tel",f)!=nothing , pyro_files )
	df_pyrohelio_files = DataFrame(:filename=>joinpath.(pyro_path,pyro_files), :start_date=>map(f->Date(match(r"neid_ljpyrohelio_chv0_(\d+).tel$",f)[1],DateFormat("yyyymmdd")),pyro_files))
end

begin
	function pick_pyrohelio_file(pyrohelio_files::DataFrame, fn::AbstractString)
		m = match(r"neidL2_(\d+T\d+).fits",fn)
		@assert m != nothing
		t_start = DateTime(m[1],DateFormat("yyyymmddTHHMMSS"))
		pick_pyrohelio_file(pyrohelio_files,t_start)
	end

	function pick_pyrohelio_file(pyrohelio_files::DataFrame, t::DateTime)
		if !issorted(pyrohelio_files, :start_date)
			sort!(pyrohelio_files, :start_date)
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
	function read_pyrohelio_file(fn::AbstractString; force::Bool = false)
                global pyrohelio_cache
		if !force && haskey(pyrohelio_cache,fn)
			df_pyrohelio_data = pyrohelio_cache[fn]
		else
			df_pyrohelio_data = CSV.read(fn, header=[:time_str,:Voltage,:FluxDensity], types=[String,Float64,Float64], DataFrame)

			datefmt = DateFormat("yyyy-mm-ddTHH:MM:SS.sss")
			df_pyrohelio_data.Time = map(tstr->DateTime(view(tstr,1:23),datefmt),collect(df_pyrohelio_data.time_str))
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
	function read_pyrohelio_file(pyrohelio_files::DataFrame, fn::AbstractString)
		m = match(r"neidL\d_(\d+T\d+).fits",fn)
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
		idx_start = searchsortedfirst(df.Time,t_start)
		if !(1 <= idx_start <= size(df,1))
			throw("get_pyrohelio_data: idx_start = " * string(idx_start))
		end
		idx_stop = searchsortedlast(df.Time,t_stop)
		if !(1 <= idx_stop <= size(df,1))
			throw("get_pyrohelio_data: idx_stop = " * string(idx_stop))
		end
		#delta_idx = searchsortedlast(view(pyrohelio_data_from_prev_file.time,idx_start:size(pyrohelio_data_from_prev_file,1)),t_stop)
		#idx_stop = idx_start+delta_idx-1
		df[idx_start:idx_stop,:]
	end

	function get_pyrohelio_data(pyrohelio_files::DataFrame, fn::AbstractString, exptime::Real)
		m = match(r"neidL\d_(\d+T\d+).fits",fn)
		@assert m != nothing
		t_start = DateTime(m[1],DateFormat("yyyymmddTHHMMSS"))
		#=
                fits = FITS(fn)
		@assert length(fits) >= 1
		hdr = read_header(fits[1])
		exptime = hdr["EXPTIME"]
		close(fits)
                =#
		exptime_sec = round(Int64,exptime)
		get_pyrohelio_data(pyrohelio_files,t_start, exptime_sec)
	end

end

function get_pyrohelio_mean_Δt(df::DataFrame; time_start::DateTime = first(df.Time) )
	@assert hasproperty(df,:Time)
	@assert hasproperty(df,:FluxDensity)
	dt = map(dt->dt.value/1000,df.Time .-time_start)
	sum(dt.*df.FluxDensity)/sum(df.FluxDensity)
end

begin
	function get_pyrohelio_summary(pyrohelio_files::DataFrame, fn::AbstractString, exptime::Real)
		m = match(r"neidL\d_(\d+T\d+).fits",fn)
		@assert m != nothing
		time_start = DateTime(m[1],DateFormat("yyyymmddTHHMMSS"))
		#=
        	fits = FITS(fn)
			@assert length(fits) >= 1
			hdr = read_header(fits[1])
			exptime = hdr["EXPTIME"]
			close(fits)
        	=#
		exptime_sec = round(Int64,exptime)
		#println(fn,": exptime = ", exptime_sec)
		try
			df = get_pyrohelio_data(pyrohelio_files, time_start, exptime_sec)
			#println(fn, ": pyrheliometer_flux = ", size(df.FluxDensity,1))
			#time_start = first(df.time)  # Comment out so use FITS file start as reference time
			mean_Δt = get_pyrohelio_mean_Δt(df, time_start=time_start)
			mean_flux = mean(df.FluxDensity)
			rms_flux = sqrt(var(df.FluxDensity,mean=mean_flux,corrected=false))
			(min_flux, max_flux) = extrema(df.FluxDensity)
			(;filename=basename(fn), time_start, exptime, mean_Δt, mean_pyroflux=mean_flux, rms_pyroflux=rms_flux, min_pyroflux=min_flux, max_pyroflux=max_flux)
		catch exepction
			(;filename=basename(fn), time_start, exptime, mean_Δt=missing, mean_pyroflux=missing, rms_pyroflux=missing,  min_pyroflux=missing, max_pyroflux=missing)
		end
	end

	function get_pyrohelio_summary(fn::String) #, exptime::Real)
		m = match(r"neidL\d_(\d+T\d+)\.fits",fn)
		@assert m != nothing
		t_start = DateTime(m[1],DateFormat("yyyymmddTHHMMSS"))
        println("# get_pyrohelio_summary on ", fn,"...") 
        flush(stdout)
		try
		println("#    Open FITS...") 
        flush(stdout)
			f = FITS(fn)
			@assert length(f) >= 1
		println("#    read_header(f[1])...") 
        flush(stdout)
			hdr = read_header(f[1])
			exptime = hdr["EXPTIME"]
		println("#    read_header(f[SOLAR])...") 
        flush(stdout)
			hdr_pyrhelio = read_header(f["SOLAR"])
			stel = hdr_pyrhelio["STEL"]
			staz = hdr_pyrhelio["STAZ"]
			#println(fn,": exptime = ", exptime)
		println("#    read Time...") 
        flush(stdout)

			t = read(f["SOLAR"],"Time")
        println("#    read Voltage...") 
        flush(stdout)
        	v = read(f["SOLAR"],"Voltage")
        println("#    read FluxDensity...") 
        flush(stdout)
        
        	fd = read(f["SOLAR"],"FluxDensity")
        	df_tmp = DataFrame(:Time=> DateTime.(t), :Voltage=>v, :FluxDensity=>fd)
			#println(fn, ": pyrheliometer_flux = ", size(df_tmp.FluxDensity,1))
			#df_tmp.Time = DateTime.(df_tmp.Time)
        	println("#    get_pyrohelio_mean_Δt...") 
            flush(stdout)
            mean_Δt = get_pyrohelio_mean_Δt(df_tmp)
			mean_pyroflux = mean(df_tmp.FluxDensity)
        	rms_pyroflux = sqrt(var(df_tmp.FluxDensity, corrected=false))
        	(min_pyroflux,max_pyroflux) = extrema(df_tmp.FluxDensity)
            println("#    close FITS...") 
            flush(stdout)
			close(f)
            println("# Finished get_pyrohelio_summary on ", fn,"...") 
            flush(stdout)
		    
			return (;filename=basename(fn), time_start=t_start, exptime, mean_Δt, mean_pyroflux, rms_pyroflux, min_pyroflux, max_pyroflux, stel, staz)
		catch exepction
            println("# Exception in get_pyrohelio_summary on ", fn,"...") 
            flush(stdout)
		
			return (;filename=basename(fn), time_start=t_start, exptime=missing, mean_Δt=missing, mean_pyroflux=missing, rms_pyroflux=missing,  min_pyroflux=missing, max_pyroflux=missing, stel=missing, staz=missing)
		end
	end

end

export make_pyrohelio_file_dataframe, pick_pyrohelio_file, read_pyrohelio_file
export get_pyrohelio_data, get_pyrohelio_summary
end # module
