#!/bin/sh

# execute in this shell 
. config

# check for CODEDIR DATADIR
if [ $CODEDIR == NULL ]; then
    echo "Error: CODEDIR is not set"
    exit 1
fi
if [ $DATADIR == NULL ]; then
    echo "Error: DATADIR is not set"
    exit 1
fi

# usage message
usage="Usage: do_stats.sh [-m model_file] "
usage="$usage [-s sim_file] [-n num_sims] [-d timestep] [-q] [model_options]"

# defaults
force=n
quiet=n

# parse argument list
while [ $# -ge 1 ]; 
  do
  
  case $1 
      in
      -m) shift; model=$1;;
      -s) shift; sim=$1;;
      -n) shift; num_sims=$1;;
      -d) shift; dt=$1;;
      -f) force="y";;
      -q) quiet="y";;
      -h) echo $usage; exit 0;;
      *) options="$options $1";;
  esac
  shift

done

# check num_sims
if [ $num_sims -gt 8999 ]; then
    echo "WARNING: num_sims > 8999, setting to 8999"
    set num_sims=8999
fi

# print message
echo "Running $num_sims simulations..."

# setting sim_options
sim_options=""

if [ $model ]; then
    sim_options="$sim_options -m $model"
fi
if [ $sim ]; then
    sim_options="$sim_options -p $sim"
fi

# set sim_command
sim_base="$CODEDIR/graph/bin/simulate"
sim_base="$sim_base $sim_options --print-stats"
sim_base="$sim_base --no-pairs --no-graph --no-degree-dist $options"
echo $sim_base

# execute 
echo -----
$sim_base
echo -----

# loop over num_sims
if [ $num_sims -gt 1 ]; then
    for ((i=101;i<(100+$num_sims);i++)); do
        $sim_base
        echo -----
    done
fi

echo 
