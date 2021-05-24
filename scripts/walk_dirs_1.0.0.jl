neid_data_path = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/"
proc_version_path = "DRPmaster/solar"
path_start = joinpath(neid_data_path,proc_version_path)
output_path_base = "/gpfs/group/ebf11/default/ebf11/neid_solar/data/outputs_1.0.0"
proj_dir = "/gpfs/group/ebf11/default/ebf11/neid_solar/code/NeidSolarScripts.jl"
script = joinpath(proj_dir,"scripts","make_manifest_solar_1.0.0.jl")

sleep_interval = 5

using CSV, DataFrames

dirs_to_process_file = joinpath(output_path_base,"dirs_to_process.csv")
if isfile(dirs_to_process_file)
   println("# Reusing ", dirs_to_process_file)
   df = CSV.read(dirs_to_process_file, DataFrame)
else
   df = DataFrame(:input_dir=>String[],:output_dir=>String[])

   for (root, dirs, files) in walkdir(path_start)
      if !occursin("level1",root)  continue end 
      for dir in dirs
              if !occursin(dir,root) continue end
              println(dir) # path to directories
              #println(joinpath(root, dir)) # path to directories
              output_dir = joinpath(output_path_base,dir)
              if !(isdir(output_dir)||islink(output_dir))
                 mkdir(output_dir)
              end
              #continuum_dir = joinpath(output_dir,"continuum")
              #if !(isdir(continuum_dir)||islink(continuum_dir))
              #   mkdir(continuum_dir)
              #end
              push!(df,(input_dir=joinpath(root,dir), output_dir=output_dir))
      end
      #=
      println("Files in $root")
      for file in files
          println(joinpath(root, file)) # path to files
      end
      =#
      # break
   end
   CSV.write(joinpath(output_path_base,"dirs_to_process.csv"),df)
end 

last_job_id = 0
function gen_job_name(; prefix::String = "auto")
  global last_job_id
  last_job_id += 1
  return prefix * string(last_job_id)
end
 
function make_pbs_scr(;proj_dir::String, script::String, input_dir::String, output_dir::String, job_name::String = gen_job_name() )
   pbs_scr = """
#!/bin/bash
#PBS -N $job_name
#PBS -l nodes=1:ppn=1
#PBS -l pmem=2000mb
#PBS -l walltime=4:00:00
###PBS -A ebf11_c_g_vc_default
###PBS -q hprc
#PBS -A cyberlamp
#PBS -l feature=rhel7
#PBS -j oe
#PBS -M ebf11@psu.edu

# Get started
echo Job started on `hostname` at `date`

# Go to the correct place
#cd \$PBS_O_WORKDIR
cd $proj_dir

# Run the job itself
~/julia1.6 --project=$proj_dir -e 'target_subdir=\"$input_dir\"; output_dir=\"$output_dir\";  include(\"$script\")'

# Finish up
echo Job Ended at `date`
"""
  return pbs_scr
end

open("submit_make_manifests.sh","w") do f_submit

for row in eachrow(df)
   # #=
    if isfile(joinpath(row.output_dir, "manifest.csv"))  
       println("# Skipping ", row.output_dir)
       continue 
    end
   # =#
    m = match(r"(\d+)$", row.input_dir)
    dir = m.captures[1] 

    pbs_scr = make_pbs_scr(proj_dir=proj_dir, script=script, input_dir=row.input_dir, output_dir=row.output_dir, job_name = "manif_"*dir)

    #println("Created script: ")
    #println(pbs_scr)
    #println("Echoing script: ")
    open("make_manifest_$dir.pbs","w") do f_pbs
       print(f_pbs, pbs_scr)
    end
    println(f_submit, "qsub make_manifest_$dir.pbs")
    println(f_submit, "sleep $sleep_interval")
    #output = readchomp(pipeline(`echo $pbs_scr`, `qsub`))
    #println("output = ")
    #println(output)
end

end # submit_make_manifests.sh

