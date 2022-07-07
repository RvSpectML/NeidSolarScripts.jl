import Pkg
Pkg.activate(mktempdir())
Pkg.add(["HTTP","EzXML"])

using Downloads
using HTTP, EzXML
pyrhelio_dir = "/gpfs/group/ebf11/default/pipeline/data/neid_solar/pyrheliometer"
pyrhelio_dir = "/storage/group/ebf11/default/pipeline/neid_solar/data/pyrheliometer"
url = "https://neid.ipac.caltech.edu/pyrheliometer.php"

try

response = HTTP.request("GET",url)
body = String(response.body)
idx = findall(r"https\:\/\/\S+\/pyroheliotel\/neid_\w+_\S{4}_\d{8}\.tel",body)
for i in eachindex(idx)
   this_url = body[idx[i]]
   m = match(r"neid_l\S+\.tel",this_url)
   if isnothing(m) continue end
   filename = joinpath(pyrhelio_dir,string(m.match))
   if filesize(filename) > 0  continue  end
   println("# Need to download ", filename)
   Downloads.download(this_url, filename)
   if !(filesize(filename) > 0)
      println("# Somehow ", filename, " is still missing/empty.") 
   end
end

catch ex
   @warn "Error: " * string(ex)
end
