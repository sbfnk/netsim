#!/bin/bash

# set script path
graph_dir=$CODEDIR/graph

gp_sim_script=$graph_dir/scripts/sim_only.gp

# check command-line args
if [[ $# != 1 ]]; then
    echo "Usage: ./makeplot.bash file_id"
    echo "       file_id is missing."
    exit 1
fi

# output directory
output_dir=$DATADIR/$1

# create gnuplot script for simulation
cat $output_dir/$1.gp > $output_dir/$1.sim.tmp.gp
sed -e s:FILE_ID:$output_dir/$1:g $gp_sim_script >> $output_dir/$1.sim.tmp.gp

# make figs
gnuplot $output_dir/$1.sim.tmp.gp > /dev/null
rm -f $output_dir/$1.sim.tmp.gp

# ps -> eps
ps2epsi sim.ps sim.eps

mv sim.eps $output_dir/$1.sim.eps

# clean
#rm -f sim.ps sim.eps
#rm -f sim.ps sim_sir.ps

exit 0
