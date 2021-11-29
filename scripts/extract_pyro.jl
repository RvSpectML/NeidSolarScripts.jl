using ArgParse
using Dates
using Base.Filesystem

function parse_commandline_extract_pyrohelio()
     s = ArgParseSettings( description = "Make pyrohelio.csv from L0 FITS files for one day.")
     @add_arg_table! s begin
         "date"
            help = "Date YYYYMMDD to extract pyroheliometer data for."
            arg_type = String
            required = true
         "--user"
            help = "NExScI Archive username"
            arg_type = String
         "--password"
            help = "NExScI Archive password"
            arg_type = String
         "--output"
            help = "Path for outputs: pyrohelio.csv"
            arg_type = String
            default = joinpath(pwd(),"pyrohelio.csv")
         "--work"
            help = "Working directory for NExScI query and L0 downloads."
            arg_type = String
         "--cookie"
            help = "Filename for NExScI cookie."
            arg_type = String
            default = "neidadmincookie.txt"
         "--query"
            help = "Filename for L0 query results."
            arg_type = String
            default = "meta_l0.csv"
         "--pyrohelio"
            help = "Path for pyroheliometer inputs: *.tel"
            arg_type = String
            #default = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/pyrohelio/202110_fromChad/"
         "--verbose"
            help = "Verbose outputs"
            action = :store_true
      end

  return parse_args(s)
end
args = parse_commandline_extract_pyrohelio()
m = match(r"(\d{4})(\d{2})(\d{2})", args["date"])
year  = parse(Int64,m[1])
month = parse(Int64,m[2])
day   = parse(Int64,m[3])
@assert 2020 <= year <= 2040
@assert 1 <= month <= 12
@assert 1 <= day <= 31

if haskey(args,"output")     output_fn          = args["output"]      end 

tmp_path = !isnothing(args["work"]) ? args["work"] : Filesystem.mktempdir()
cookie_fn          = args["cookie"]   
query_fn           = args["query"]   
pyrohelio_dir      = args["pyrohelio"] 
verbose = args["verbose"]

using FITSIO
using CSV, DataFrames
using Statistics
using NeidArchive

cookie_nexsci = joinpath(pwd(),cookie_fn) # "./neidadmincookie.txt"
query_result_file = joinpath(pwd(),query_fn) # "./criteria1.csv"
println("tmp_path = >",tmp_path,"<.")
#cookie_nexsci = joinpath(tmp_path,cookie_fn) # "./neidadmincookie.txt"
#query_result_file = joinpath(tmp_path,query_fn) # "./criteria1.csv"

if !(filesize(cookie_nexsci)>0)
   if !isnothing(args["user"]) && !isnothing(args["password"])
      NeidArchive.login(userid=args["user"], password=args["password"], cookiepath=cookie_nexsci)
   else
      NeidArchive.login(cookiepath=cookie_nexsci)
   end
end

param = Dict{String,String}()
 param["datalevel"] = "solarl0"
 param["piname"] = "Mahadevan"
 param["object"] = "Sun"
 param["datetime"] = NeidArchive.datetime_one_day_solar(Date(year,month,day))

NeidArchive.query(param, cookiepath=cookie_nexsci, outpath=query_result_file)# , outdir=tmp_path)
num_lines = countlines(query_result_file) - 1
println("# Query resulted in file with ", num_lines, " entries.")

start_dir = pwd()
#cd(tmp_path)
#NeidArchive.download(query_result_file, param["datalevel"], cookiepath=cookie_nexsci, start_row=1, end_row=3)
NeidArchive.download(query_result_file, "l0", cookiepath=cookie_nexsci, outdir=tmp_path) # , start_row=1, end_row=3)

touch(output_fn)
exampleFileIOStream =  open(output_fn,"w")
#cd(start_dir)
close(exampleFileIOStream)
exit(0)

fnin = ""
for (i,line) in enumerate(lines)
   m = match(r"(neidL0_\d{4}\d{2}\d{2}T\d+\.fits)", line)
   if (m==nothing) continue end 
   fnin = string(m[1])
   fninfull = joinpath(tmp_path,fnin)
   if !(filesize(fninfull) > 0)
      m = match(r"neidL0_(\d{4})(\d{2})(\d{2})T(\d+)\.fits", fnin)
      if (m==nothing) continue end 
      println("# Downloading $fnin.\n")
      try
        yyyy = m[1]
        mm = m[2]
        dd = m[3]
        tttttt = m[4]
        run(`wget --load-cookies neid_cookie.txt -O $fninfull "https://neid.ipac.caltech.edu/get_file.php?filehand=raw/$yyyy/$mm/$dd/$fnin&solar"`)
      catch ex
         println("# Failed to download $fninfull.\n")
         continue
      end
    end
    fnout = replace(fnin,"neidL0_"=>"neid_pyro_")
    fnout = replace(fnout,"fits"=>"csv")
    df = DataFrame
    if  filesize(fnout) >0  
        df = CSV.read(fnout,DataFrame)
    else
    try
       f = FITS(fninfull)
       t = read(f[21],"Time")
       v = read(f[21],"Voltage")
       fd = read(f[21],"FluxDensity")
       df = DataFrame(:Time=>t, :Voltage=>v, :FluxDensity=>fd)
       CSV.write(fnout, df)
     catch ex
        println("# Failed to write $fnout.\n")
         continue
     end
     end
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
     rm(fninfull)
end
close(exampleFileIOStream)

