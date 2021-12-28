using Logging
using ArgParse

@info "# Parsing arguments..."
 function parse_commandline()
     s = ArgParseSettings( description = "Combined NEID solar daily reports into one csv file.")
     add_arg_group!(s, "Files", :argg_files)
     @add_arg_table! s begin
         "input_path"
             help = "Path for input files"
             arg_type = String
             default = pwd()
         "output"
             help = "Filename for output (csv)"
             arg_type = String
             default = "summary.csv"
         "--input_filename"
             help = "Filename for daily reports (toml)"
             arg_type = String
             default = "daily_summary.toml"
         "--exclude_filename"
             help = "Filename with dates to exclude (csv)"
             arg_type = String
         "--overwrite"
            help = "Specify it's ok to overwrite the output file."
            #default = true
            action = :store_true
         "--log_info"
            help = "Print info-level messages."
            action = :store_true
         "--log_debug"
            help = "Print debug-level messages."
            action = :store_true
      end

     return parse_args(s)
 end
 args = parse_commandline()
 if args["log_debug"]
    global_logger(ConsoleLogger(stderr, Logging.Debug))
 elseif args["log_info"]
    global_logger(ConsoleLogger(stderr, Logging.Info))
 end

 input_fn = args["input_filename"]
 input_path = args["input_path"]
 output_fn = joinpath(input_path,args["output"])

if !(args["overwrite"] || !isfile(output_fn) || (filesize(output_fn)==0))
   @error "Can't overwrite " output_filename=output_fn
   exit(1)
else
   #touch(output_fn)  # Should we create empty file as a lock?
end

@info "# Loading packages"
using Dates, Markdown, TOML
using CSV, DataFrames, Query, Glob

@info "# Finding files named $input_fn"
files1 = glob([r"\d{2}",args["input_filename"]],args["input_path"])
files2 = glob([r"\d{2}",r"\d{2}",args["input_filename"]],args["input_path"])
files3 = glob([r"\d{4}",r"\d{2}",r"\d{2}",args["input_filename"]],args["input_path"])
files = vcat(files1,files2,files3)
@assert length(files) >= 1

if !isnothing(args["exclude_filename"]) && isfile(args["exclude_filename"]) && (filesize(args["exclude_filename"])>0)
   @info "# Reading days to exclude."
   df_exclude = CSV.read(args["exclude_filename"],DataFrame)
   @info "# Found $(size(df_exclude,1)) dates to exclude"
else
   df_exclude = DataFrame()
end

#=
function dates_to_strings!(d::AbstractDict)
   for (k,v) in pairs(d)
      if typeof(v) <: Date
         d[k] = string(v)
      elseif typeof(v) <: DateTime
         d[k] = string(v)
      elseif typeof(v) <: AbstractDict
         dates_to_strings!(d[k])
      end
   end
   return d
end 
=#

function _flatten_dict!(d_out::D1, d_in::D2; prestr::String = "") where { D1<:AbstractDict,  D2<:AbstractDict } 
   for (k1,v1) in pairs(d_in)
       if typeof(v1) <: AbstractDict
          _flatten_dict!(d_out,v1, prestr = prestr * k1 * "." ) 
       else
          d_out[prestr * k1 ] = v1
       end
   end
   return d_out
end

function flatten_dict(d::D) where { D<:AbstractDict } 
   d_out =  D()   
   _flatten_dict!(d_out, d)
end

function array_of_dict_to_dataframe( aod::AOD; first_cols::AbstractVector{String}=String[], exclude_cols::AbstractVector{String}=String[] ) where { D<:AbstractDict, AOD<:AbstractVector{D} }
  df = DataFrame()
  for k in first_cols
      if !haskey(first(aod), k) || in(k, exclude_cols)
         continue 
      end 
      df[!,k] = map(i->aod[i][k],1:size(aod,1)) 
  end
  for k in sort!(collect(keys(first(aod))))
      if in(k, first_cols) || in(k, exclude_cols) 
         continue
      end
      df[!,k] = map(i->aod[i][k],1:size(aod,1)) 
  end
  return df
end

@info "# Parsing daily toml files."  num_files=length(files)
daily = Vector{Dict{String,Any}}(undef, size(files,1) )
println("files = ", files)
flush(stdout)
flush(stderr)
j = 0
parser = TOML.Parser()
for file in files
    if filesize(file) >0 
       @info "# Processing $file"
       d = TOML.parsefile(parser,file) 
       if haskey(d,"obs_date") && haskey(d["obs_date"],"string") && (size(df_exclude,1)>=1) && ( in(d["obs_date"]["string"], df_exclude.date_to_exclude) )
          continue
       end
       #dates_to_strings!(d)
       global j += 1
       daily[j] = d
    end
end
num_days_with_usable_obs = j
resize!(daily,num_days_with_usable_obs)

@info "# Making dataframe"
daily_flat = flatten_dict.(daily)
first_cols = ["obs_date.string", "obs_date.mean_bjd","num_rvs.good"]
cols_to_exclude = ["report.hostname", "title", "report.input_md5", "report.input_file"]
df =  array_of_dict_to_dataframe( daily_flat, first_cols=first_cols, exclude_cols=cols_to_exclude )

@info "# Sorting dataframe"
df_sorted = sort!(df,Symbol("obs_date.mean_bjd") ) 

#=
df.obs_date = map(day->day["obs_date"]["string"], daily)
df.obs_date = map(day->day["obs_date"]["string"], daily)
df.mean_bjd = map(day->day["obs_date"]["mean_bjd"], daily)
df.num_rvs_usable = map(day->day["num_rvs"]["usable"], daily)
df.num_rvs_good = map(day->day["num_rvs"]["good"], daily)
df.mean_rv_drp = map(day->day["rv"]["drp"]["mean_rv"], daily)
df.median_rv_drp = map(day->day["rv"]["drp"]["mean_rv"], daily)
df.median_σ_rv = map(day->day["rv"]["drp"]["median_σ_rv"], daily)
df.rms_rv_drp = map(day->day["rv"]["drp"]["rms_rvs"], daily)
df.winsor_mean_rv_drp = map(day->day["rv"]["drp"]["winsor_mean_rv"], daily)
df.winsor_rms_rv_drp = map(day->day["rv"]["drp"]["winsor_rms_rv"], daily)
@info "# Sorting dataframe"
df_sorted = df |> @orderby( _.mean_bjd ) |> DataFrame
=#

@info "# Writing CSV file"
CSV.write(output_fn,df_sorted )


