using PyCall
using Conda
Conda.add("astropy")
#Conda.add("astropy.coordinates")
#Conda.add("astropy.time")
#Conda.add("astroquery.jplhorizons")
Conda.add("astroquery")
Conda.add("barycorrpy")
#Conda.add("scipy.spatial.transform")
Conda.add("scipy")

using Conda
pipcmd = joinpath(Conda.PYTHONDIR,"pip")
run(`$pipcmd install pyneid`)

