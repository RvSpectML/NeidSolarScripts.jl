#!/bin/bash
export JULIA_NUM_THREADS=4
julia -e 'idx_day_to_use=1; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=2; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=3; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=4; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=5; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=6; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=7; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=8; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=9; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=10; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=11; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=12; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=13; include("examples/calc_ccfs.jl")';
julia -e 'idx_day_to_use=14; include("examples/calc_ccfs.jl")';
