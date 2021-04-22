ENV["PYTHON"] = ""
using PyCall
using Conda

pipcmd = joinpath(Conda.PYTHONDIR,"pip")
#Conda.add("astroquery.jplhorizons")
#Conda.add("astroquery")
run(`$pipcmd install --pre astroquery`)

Conda.add("astropy")
#Conda.add("astropy.coordinates")
#Conda.add("astropy.time")
Conda.add("barycorrpy")
#Conda.add("scipy.spatial.transform")
Conda.add("scipy")

#run(`$pipcmd install pyneid`)

