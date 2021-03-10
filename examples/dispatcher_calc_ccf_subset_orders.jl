#!/bin/bash
export JULIA_NUM_THREADS=4

#good_orders = 55 .+ collect(1:14);
 #good_orders = 55 .+ collect(17:24)
 #good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41))
# julia -e 'idx_day_to_use=1; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
# julia -e 'idx_day_to_use=2; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
# julia -e 'idx_day_to_use=3; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
# julia -e 'idx_day_to_use=4; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
# julia -e 'idx_day_to_use=5; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
# julia -e 'idx_day_to_use=6; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
# julia -e 'idx_day_to_use=7; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
# julia -e 'idx_day_to_use=8; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
# julia -e 'idx_day_to_use=9; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
# julia -e 'idx_day_to_use=10; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
# julia -e 'idx_day_to_use=11; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
# julia -e 'idx_day_to_use=12; good_orders = 55 .+ collect(1:14); include("examples/calc_ccfs_subset_orders.jl")';
#julia -e 'idx_day_to_use=13; include("examples/calc_ccfs.jl")';
#julia -e 'idx_day_to_use=14; include("examples/calc_ccfs.jl")';

#good_orders = 55 .+ collect(17:24);
 julia -e 'idx_day_to_use=1; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=1; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';

 julia -e 'idx_day_to_use=2; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=2; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';

 julia -e 'idx_day_to_use=3; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=3; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=4; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=4; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';

 julia -e 'idx_day_to_use=5; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=5; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';

 julia -e 'idx_day_to_use=6; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=6; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';

  julia -e 'idx_day_to_use=7; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=8; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=9; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=10; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=11; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=12; good_orders = 55 .+ collect(17:24); include("examples/calc_ccfs_subset_orders.jl")';
#julia -e 'idx_day_to_use=13; include("examples/calc_ccfs.jl")';
#julia -e 'idx_day_to_use=14; include("examples/calc_ccfs.jl")';

#good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));
 julia -e 'idx_day_to_use=7; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=8; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=9; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=10; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=11; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';
 julia -e 'idx_day_to_use=12; good_orders = 55 .+ vcat(collect(27:29), collect(37:38),collect(40:41));  include("examples/calc_ccfs_subset_orders.jl")';
#julia -e 'idx_day_to_use=13; include("examples/calc_ccfs.jl")';
#julia -e 'idx_day_to_use=14; include("examples/calc_ccfs.jl")';
