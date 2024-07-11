#project_dir = "/storage/work/ebf11/LineByLineFits"
#import Pkg
#Pkg.activate(project_dir)

@info ARGS

using HDF5, JLD2, CSV, FITSIO, FileIO
using DataFrames, Dates, InlineStrings, OrderedCollections
using StatsBase, Statistics, ProgressMeter

root_dir = "/storage/group/ebf11/default/pipeline/neid_solar"
output_dir = length(ARGS) >= 1 ? ARGS[1] : "data/v1.3/outputs/ebf11/Espresso_cont3" # "data/v1.3/outputs/ebf11/linefinder_v1.3b/")
output_dir = joinpath(root_dir,output_dir)
@info output_dir
@assert isdir(output_dir)

daily_filename = length(ARGS) >= 2 ? ARGS[2] : "daily_ccfs_espressoG2_norm=cont&mask=3.jld2"
#daily_filename = "daily_ccfs_bf=3_depth=0.02_q=95_fconv=95_telluric=10000_norm=cont&mask=3&tophat.jld2"
@info daily_filename
@assert isfile(joinpath(output_dir,"2021","01","01",daily_filename))
mean_daily_filename = "mean_" * daily_filename 
if filesize(mean_daily_filename)>0 
  @info "Desired output file already exists: ", mean_daily_filename
  exit(0)
end

function find_file_in_subdirs(dir::AbstractString; file_to_find::AbstractString)
    filelist = String[]
    for (root,dirs,files) in walkdir(dir)
        for (i,file) in enumerate(files)
            if occursin(file_to_find,file)
                push!(filelist, joinpath(root, file))
            end
        end
    end
    return filelist
end


function path_to_date(fn::AbstractString)
    m = match(r"\/(\d{4})\/(\d{2})\/(\d{2})\/",fn)
    Date(parse(Int,m[1]),parse(Int,m[2]),parse(Int,m[3]))
end


function get_num_clean_obs(fn::AbstractString)
   try 
        sum(load(fn,"clean_obs_mask"))
    catch
        @warn "Failed to read clean_obs_mask from ", fn
        return 0
    end
end

#df_days_exclude = DataFrame(:date_to_exclude=>Date[]);
df_days_exclude = CSV.read(joinpath(root_dir,"data","days_to_exclude.csv"),DataFrame)

files_to_process = find_file_in_subdirs(joinpath(output_dir), file_to_find=daily_filename)
good_files_to_process = filter(fn->filesize(fn)>0, files_to_process)
good_files_to_process = filter(fn->!occursin(mean_daily_filename,fn), good_files_to_process)
@info "Finding days with >= 100 good observations ", size(good_files_to_process)
good_files_to_process = filter(fn->get_num_clean_obs(fn)>=100, good_files_to_process)
@info "Skipping days to exclude ", size(good_files_to_process)
good_files_to_process = filter(fn->path_to_date(fn) ∉ df_days_exclude.date_to_exclude ,good_files_to_process)
@info "Skipping days to exclude ", size(good_files_to_process)
date_strs_processed = replace.(map(fn->dirname(fn)[end-9:end],good_files_to_process),"/"=>"-")
length(good_files_to_process)

df_by_date = DataFrame(:date=>Date.(date_strs_processed),:Filename=>good_files_to_process);

@info "Loading manifests"
df_by_date.manifests = map(fn->load(fn,"manifest")[!,[:Filename,:time_start,:drpextsnr,:drp_ccfjdmod,:drp_dvrmsmod,:Δv_diff_ext,:Δfwhm²,:airmass,:expmeter_mean,:expmeter_rms,:mean_pyroflux,:rms_pyroflux,:CaIIHK,:σ_CaIIHK]],df_by_date.Filename);
@info "Done"

function compute_daily_order_ccfs(df_day::AbstractDataFrame; num_orders = 122, num_vels = 1604)
    order_ccfs_day = zeros(num_vels,num_orders)
    num_obs = 0
    for fn in df_day.Filename
        try
            f = FITS(fn)
            order_ccfs = read(f["CCFS"])
            order_ccfs += order_ccfs
            num_obs += 1
        catch
            println("# Skipping fn")
        end
    end
    order_ccfs /= num_obs
end

function compute_daily_mean_ccf(fn::AbstractString)
    #v_grid_first = load(df_by_date.Filename[1],"v_grid")
    #orders_to_use_first = load(df_by_date.Filename[1],"orders_to_use")
    clean_obs_mask =  load(fn,"clean_obs_mask")
    num_obs_clean = sum(clean_obs_mask)
    order_ccfs_first = load(fn,"order_ccfs")
    order_ccf_vars_first = load(fn,"order_ccf_vars")
    #size(order_ccfs_first), size(clean_obs_mask), size(v_grid_first), orders_to_use_first
    mean_order_ccf = mean(view(order_ccfs_first,:,:,clean_obs_mask),dims=3)
    mean_order_ccf_var = mean(view(order_ccf_vars_first,:,:,clean_obs_mask),dims=3)./sqrt(num_obs_clean)
    return (; mean_order_ccf, mean_order_ccf_var, num_obs_clean)
end

#(ccf_flux_tmp, ccf_var_tmp, num_obs_clean_tmp) = compute_daily_mean_ccf(df_by_date.Filename[2])

function compute_daily_order_ccfs(fn_list::Vector{String}) # ; num_orders = 122, num_vels = 1604)
    num_days = length(fn_list)
    (ccf_flux_first, ccf_var_first) = compute_daily_mean_ccf(first(fn_list))
    num_vels = size(ccf_flux_first,1)
    num_orders = size(ccf_flux_first,2)
    daily_order_ccfs = zeros(num_vels,num_orders,num_days)
    daily_order_ccf_vars = zeros(num_vels,num_orders,num_days)
    num_obs_clean = zeros(num_days)
    @showprogress for (i,fn) in enumerate(fn_list)
        #if mod(i,5)==1    println("# Reading " * fn)    end
        #try
            (ccf_flux_tmp, ccf_var_tmp, num_obs_clean_tmp) = compute_daily_mean_ccf(fn)
            daily_order_ccfs[:,:,i] .= ccf_flux_tmp
            daily_order_ccf_vars[:,:,i] .= ccf_var_tmp
            num_obs_clean[i] = num_obs_clean_tmp
        #=catch
            println("# Skipping " *fn)
            continue
        end =#
    end
    return (;daily_order_ccfs, daily_order_ccf_vars, num_obs_clean)
end

@info "Computing daily order CCFs"
(daily_order_ccfs, daily_order_ccf_vars, num_obs_clean) = compute_daily_order_ccfs(df_by_date.Filename)

@info "Reading v_grid"
v_grid_first = load(df_by_date.Filename[1],"v_grid")
@info "Reading Orders to use"
orders_to_use_first = load(df_by_date.Filename[1],"orders_to_use")

@info "Saving results"
jldsave(joinpath(output_dir,mean_daily_filename); daily_order_ccfs, daily_order_ccf_vars, v_grid=v_grid_first, orders_to_use= orders_to_use_first, dates=df_by_date.date, num_obs_clean, Filename = df_by_date.Filename)
    #manifest=df_by_date.manifests, dates=df_by_date.date, num_obs_day = map(df->size(df,1),df_by_date.manifests) )

@info "Wrote " mean_daily_filename
