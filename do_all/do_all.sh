#!/bin/sh

. $CODEDIR/do_all/config

usage="Usage: do_all.sh file_id [-m model_file] [-o ode_file]"
usage="$usage [-s sim_file] [-n num_sims] [-d timestep] [model_options]"

file_id=$1
output_dir="$DATADIR/$file_id"

force=n

if [ -z $file_id ]; then
  echo "Error: file_id not set"
  echo "$usage"
  exit 1
fi

shift

# parse argument list

while [ $# -ge 1 ]; 
do
  
  case $1 
  in
      -m) shift; model=$1;;
      -o) shift; ode=$1;;
      -s) shift; sim=$1;;
      -n) shift; num_sims=$1;;
      -d) shift; dt=$1;;
      -f) force="y";;
      -h) echo $usage; exit 0;;
      *) options="$options $1";;
  esac
  shift

done

if [ -d $output_dir ]; then
    if [ $force == "y" ]; then
	echo "Overwriting contents of $output_dir"
        if [ -d $output_dir/images ]; then
  	    rm -f $output_dir/images/*
        else
            mkdir $output_dir/images
        fi
    else
	echo "Error: $output_dir exists"
	echo "Use -f to override"
	exit 1
    fi
else
    mkdir $output_dir
    mkdir $output_dir/images
fi

if [ $num_sims -gt 8999 ]; then
  echo "WARNING: num_sims > 8999, setting to 8999"
  set num_sims=8999
fi

echo "Running $num_sims simulations..."

sim_options=""

if [ $model ]; then
    sim_options="$sim_options -m $model"
fi
if [ $sim ]; then
    sim_options="$sim_options -p $sim"
fi

sim_base="$CODEDIR/graph/simulate"
sim_base="$sim_base --graph-dir $output_dir/images $sim_options"
sim_base="$sim_base $options"
sim_command="$sim_base --write-file $output_dir/$file_id""100"

$sim_command

if [ $num_sims -gt 1 ]; then
    for ((i=101;i<(100+$num_sims);i++)); do
	sim_command="$sim_base --write-file $output_dir/$file_id$i"
        $sim_command
    done
fi

echo "Averaging runs..."

avg_command="$CODEDIR/graph/average_runs -d $dt -o $output_dir/$file_id"
avg_command="$avg_command $output_dir/$file_id???.sim.dat"

$avg_command

rm $output_dir/$file_id???.sim.dat

echo "Running ode..."

if [ $model ]; then
  ode_options="$ode_options -m $model"
fi
if [ $ode ]; then
  ode_options="$ode_options -p $ode"
fi

ode_command="$CODEDIR/ode/ode_solve --file_id $output_dir/$file_id"
ode_command="$ode_command $ode_options $options"
$ode_command

echo "Making plots"
plot_command="$CODEDIR/do_all/makeplot.bash $file_id"
$plot_command
