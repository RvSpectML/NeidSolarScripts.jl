### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ fae2ba1d-02ab-4cee-8c83-abda1e6a098c
begin
	import Pkg
	#=Pkg.activate(mktempdir())
	Pkg.add([
			Pkg.PackageSpec(name="CSV"), 
			Pkg.PackageSpec(name="DataFrames")#,
			#Pkg.PackageSpec(name="EchelleCCFs")
			])
	Pkg.develop(name="EchelleCCFs")
	=#using CSV, DataFrames
	using EchelleCCFs 
end

# ╔═╡ 0cbdbb36-4874-4862-823f-6bf0fbd4e60c
begin
	neid_data_path = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/"
	proc_version_path = "DRPmaster/solar"
	path_start = joinpath(neid_data_path,proc_version_path)
	output_path_base = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/outputs_1.0.0"
	proj_dir = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl"
	script = joinpath(proj_dir,"examples","calc_order_ccfs_using_continuum_1.0.0.jl")
end

# ╔═╡ a62db4ef-c321-4dfe-b167-19c29ce45484
html"""<style>
main {
    max-width: 1000px;
}
"""

# ╔═╡ bd273401-b32e-44a6-a490-92575ed05ae9
begin 
	df = DataFrame(:manifest_filename=>String[], :output_filename=>String[],:num_lines=>Int64[])
	for (root, dirs, files) in walkdir(output_path_base)
		for file in files
			if !occursin(r"manifest\.csv",file)  continue end 
		  	#println("root = ", root)
		  	#println("dirs = ", dirs)
			#println("files = ", files)
			manifest_fn = joinpath(root,file)
			output_fn = joinpath(root,"daily_ccfs.jld2")
			nlines = countlines(manifest_fn)
			push!(df, Dict(:manifest_filename=>manifest_fn, :output_filename=>output_fn, :num_lines=>nlines))
    	end
	end
	df
end

# ╔═╡ d48ad3e0-14db-47b4-84db-4db10e1f0114
begin
	num_threads = 1
	min_order = 56
	max_order = 108
	#line_list_filename = joinpath(ENV["JULIA_DEPOT_PATH"],"dev/EchelleCCFs/data/masks","espresso+neid_mask_97_to_108.mas")
        #line_list_filename = joinpath(pkgdir(EchelleCCFs),"data/masks","espresso+neid_mask_97_to_108.mas")
        line_list_filename1 = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl/scripts/linelist_20210208.csv"
        line_list_filename2 = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl/data/VALD_species=all_depthcutoff=0.05_overlapcutoff=1e-05_allowBlends=0_badLineFilter=20210115_solar_rejectTelluricSlope=0.0_waves=Reiners_depths=original_nbin=1_binparam=depth_n=0.csv"
        #line_list_filename3 = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl/data/VALD_species=all_depthcutoff=0.05_overlapcutoff=1e-05_allowBlends=0_badLineFilter=ESPRESSOG2_rejectTelluricSlope=0.0_waves=Reiners_depths=original_nbin=1_binparam=depth_n=0.csv"
        #line_list_filename4 = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl/data/VALD_species=all_depthcutoff=0.05_overlapcutoff=2e-05_allowBlends=0_badLineFilter=20210115_solar_rejectTelluricSlope=0.0_waves=Reiners_depths=original_nbin=1_binparam=depth_n=0.csv"
        #line_list_filename5 = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl/data/VALD_species=all_depthcutoff=0.05_overlapcutoff=2e-05_allowBlends=0_badLineFilter=ESPRESSOG2_rejectTelluricSlope=0.0_waves=Reiners_depths=original_nbin=1_binparam=depth_n=0.csv"
        #line_list_filename6 = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl/data/VALD_species=all_depthcutoff=0.05_overlapcutoff=3e-05_allowBlends=0_badLineFilter=20210115_solar_rejectTelluricSlope=0.0_waves=Reiners_depths=original_nbin=1_binparam=depth_n=0.csv"
        #line_list_filename7 = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl/data/VALD_species=all_depthcutoff=0.05_overlapcutoff=3e-05_allowBlends=0_badLineFilter=ESPRESSOG2_rejectTelluricSlope=0.0_waves=Reiners_depths=original_nbin=1_binparam=depth_n=0.csv"
        #line_list_filename8 = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl/data/VALD_species=all_depthcutoff=0.05_overlapcutoff=6e-05_allowBlends=0_badLineFilter=20210115_solar_rejectTelluricSlope=0.0_waves=Reiners_depths=original_nbin=1_binparam=depth_n=0.csv"
        #line_list_filename9 = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl/data/VALD_species=all_depthcutoff=0.05_overlapcutoff=6e-05_allowBlends=0_badLineFilter=ESPRESSOG2_rejectTelluricSlope=0.0_waves=Reiners_depths=original_nbin=1_binparam=depth_n=0.csv"
        line_list_filename_list = [line_list_filename1, line_list_filename2 ]
        #line_list_filename_list = [line_list_filename1, line_list_filename2, line_list_filename3, line_list_filename4, line_list_filename5, line_list_filename6, line_list_filename7, line_list_filename8 ]
	sed_filename = joinpath(proj_dir,"data", "neidMaster_HR_SmoothLampSED_20210101.fits")	
        #anchors_filename = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl/scripts/anchors_20210208.jld2"
        anchors_filename = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl/scripts/anchors_57_58.jld2"
end

# ╔═╡ 8366a29a-1857-410f-aea3-87049e6248f6
pkgdir(EchelleCCFs)

# ╔═╡ e08192ce-9f0d-4e6f-a4aa-864f093dfcdd
#df[1,:manifest_filename]
df[1,:output_filename]

# ╔═╡ 73d682e9-cfa2-45b5-ba9e-b7cbf157b31d
begin
	last_job_id = 0
	function gen_job_name(; prefix::String = "auto")
	  global last_job_id
	  last_job_id += 1
	  return prefix * string(last_job_id)
	end
end

# ╔═╡ ed421e9d-2d2a-4af6-b832-591f38e092ad
function make_pbs_scr(;proj_dir::String, script::String, input_fn::String, output_fn::String, job_name::String = gen_job_name(),  min_row::Integer = 1, max_row::Integer = size(df,1)  )
	pbs_hdr_str = """
#!/bin/bash
#PBS -N $job_name
#PBS -l nodes=1:ppn=$num_threads
#PBS -l pmem=4000mb
#PBS -l walltime=2:00:00
###PBS -A ebf11_c_g_vc_default
###PBS -q hprc
#PBS -A cyberlamp
#PBS -l feature=rhel7
#PBS -j oe
#PBS -M ebf11@psu.edu

# Get started
echo Job started on `hostname` at `date`

# Go to the correct place
cd \$PBS_O_WORKDIR
#cd $proj_dir
"""

        pbs_str = pbs_hdr_str
        for (i,row) in enumerate(eachrow(df))
                if !(min_row<=i<=max_row) continue end
                input_fn = df[i,"manifest_filename"]
                output_fn = df[i,"output_filename"]

                for (j,line_list_filename) in enumerate(line_list_filename_list)

                   output2_fn = replace(output_fn,"daily_ccfs"=>"daily_ccfs_" * string(j) )
                   pbs_cmd_str = """
# Run the job itself
~/julia --project=$proj_dir -t $num_threads $script $input_fn $output2_fn --line_list_filename $line_list_filename  --sed_filename $sed_filename  --anchors_filename $anchors_filename  --orders_to_use=$min_order $max_order --range_no_mask_change 6.0 --apply_continuum_normalization  --variable_mask_scale --overwrite
"""
                  pbs_str = pbs_str * pbs_cmd_str
                end
        end

pbs_ftr_str = """
# Finish up
echo Job Ended at `date`
"""
	pbs_str = pbs_str * pbs_ftr_str
  return pbs_str
end


# ╔═╡ b5769171-0b0d-4e48-ae4e-be8c4995b4f2
"output_filename" ∈ names(df)

# ╔═╡ 8a7ff512-39aa-4b76-b3fd-f02ac410c85b

function make_pbs_multi_scr(;proj_dir::String, script::String, df::DataFrame, job_name::String = "auto", num_threads::Integer=1, pmem::String ="8000mb", walltime::String ="23:59:00", min_row::Integer = 1, max_row::Integer = size(df,1) ) 
	@assert "manifest_filename" ∈ names(df)
	@assert "output_filename" ∈ names(df)
    @assert size(df,1) >=1
	
	pbs_hdr_str = """
#!/bin/bash
#PBS -N $job_name
#PBS -l nodes=1:ppn=$num_threads
#PBS -l pmem=$pmem
#PBS -l walltime=$walltime
###PBS -A ebf11_c_g_vc_default
###PBS -q hprc
#PBS -A cyberlamp
#PBS -l feature=rhel7
#PBS -j oe
#PBS -M ebf11@psu.edu

# Get started
echo Job started on `hostname` at `date`

# Go to the correct place
cd \$PBS_O_WORKDIR
#cd $proj_dir
"""
 	pbs_str = pbs_hdr_str
	for (i,row) in enumerate(eachrow(df))
		if !(min_row<=i<=max_row) continue end
		input_fn = df[i,"manifest_filename"]
		output_fn = df[i,"output_filename"]

                for (j,line_list_filename) in enumerate(line_list_filename_list)

                   output_fn2 = replace(output_fn,"daily_ccfs"=>"daily_ccfs_" * string(j) )
   		   pbs_cmd_str = """
# Run the job itself
~/julia --project=$proj_dir -t $num_threads $script $input_fn $output_fn2 --line_list_filename $line_list_filename  --sed_filename $sed_filename  --anchors_filename $anchors_filename  --orders_to_use=$min_order $max_order --range_no_mask_change 6.0  --apply_continuum_normalization  --variable_mask_scale  --overwrite

"""
		   pbs_str = pbs_str * pbs_cmd_str
                end

	end
pbs_ftr_str = """
	# Finish up
echo Job Ended at `date`
"""

	pbs_str = pbs_str * pbs_ftr_str
  return pbs_str
end


# ╔═╡ 962d7e9e-4fdd-4080-82bf-f0175002ca19

open("submit_calc_ccfs.sh","w") do f_submit

for row in eachrow(df)
	if row.num_lines <= 1 continue end
   #=
    if isfile(joinpath(row.output_dir, "manifest.csv"))  
       println("# Skipping ", row.output_dir)
       continue 
    end
   =#
    m = match(r"(\d+)$", dirname(row.manifest_filename))
    dir = m.captures[1] 

    pbs_scr = make_pbs_scr(proj_dir=proj_dir, script=script, input_fn=row.manifest_filename, output_fn=row.output_filename, job_name = "ccfs_"*dir)

    #println("Created script: ")
    #println(pbs_scr)
    #println("Echoing script: ")
    open("ccfs_$dir.pbs","w") do f_pbs
       print(f_pbs, pbs_scr)
    end
    println(f_submit, "qsub ccfs_$dir.pbs")
    #println(f_submit, "sleep $sleep_interval")
    #output = readchomp(pipeline(`echo $pbs_scr`, `qsub`))
    #println("output = ")
    #println(output)
end

end # submit_calc_ccfs.sh



# ╔═╡ 487c2c86-e56f-4823-be3c-05c076dbbe89

open("submit_calc_ccfs_multi.sh","w") do f_submit
files_per_jobs = 8
for i in reverse(1:files_per_jobs:size(df,1))
	#if row.num_lines <= 1 continue end
   #=
    if isfile(joinpath(row.output_dir, "manifest.csv"))  
       println("# Skipping ", row.output_dir)
       continue 
    end
   =#
    #m = match(r"(\d+)$", dirname(row.manifest_filename))
    #dir = m.captures[1] 

    pbs_scr = make_pbs_multi_scr(proj_dir=proj_dir, script=script, df=df, min_row=i, max_row=min(i+files_per_jobs-1,size(df,1)), job_name = "ccfs_"*string(i), num_threads=num_threads)

    #println("Created script: ")
    #println(pbs_scr)
    #println("Echoing script: ")
    open("ccfs_" *string(i) * ".pbs","w") do f_pbs
       print(f_pbs, pbs_scr)
    end
    println(f_submit, "qsub ccfs_" *string(i) *".pbs")
    #println(f_submit, "sleep $sleep_interval")
    #output = readchomp(pipeline(`echo $pbs_scr`, `qsub`))
    #println("output = ")
    #println(output)
end

end # submit_calc_ccfs.sh|


# ╔═╡ Cell order:
# ╠═fae2ba1d-02ab-4cee-8c83-abda1e6a098c
# ╠═0cbdbb36-4874-4862-823f-6bf0fbd4e60c
# ╠═a62db4ef-c321-4dfe-b167-19c29ce45484
# ╠═bd273401-b32e-44a6-a490-92575ed05ae9
# ╠═d48ad3e0-14db-47b4-84db-4db10e1f0114
# ╠═8366a29a-1857-410f-aea3-87049e6248f6
# ╠═e08192ce-9f0d-4e6f-a4aa-864f093dfcdd
# ╠═73d682e9-cfa2-45b5-ba9e-b7cbf157b31d
# ╠═ed421e9d-2d2a-4af6-b832-591f38e092ad
# ╠═b5769171-0b0d-4e48-ae4e-be8c4995b4f2
# ╠═8a7ff512-39aa-4b76-b3fd-f02ac410c85b
# ╠═962d7e9e-4fdd-4080-82bf-f0175002ca19
# ╠═487c2c86-e56f-4823-be3c-05c076dbbe89
