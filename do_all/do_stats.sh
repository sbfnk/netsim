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
usage="$usage [-s sim_file] [-d timestep] [-q] [model_options]"

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
      -d) shift; dt=$1;;
      -f) force="y";;
      -q) quiet="y";;
      -h) echo $usage; exit 0;;
      *) options="$options $1";;
  esac
  shift

done

# setting sim_options
sim_options=""

if [ $model ]; then
    sim_options="$sim_options -m $model"
fi
if [ $sim ]; then
    sim_options="$sim_options -p $sim"
fi

# set sim_command
sim_base="$CODEDIR/sim/bin/simulate"
sim_base="$sim_base $sim_options --print-stats"
sim_base="$sim_base $options"
echo $sim_base

# execute 
$sim_base

echo 
