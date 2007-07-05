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
usage="Usage: do_sim.sh file_id [-m model_file] "
usage="$usage [-s sim_file] [-d timestep] [-q] [model_options]"

# set file_id
file_id=$1

# set output dir
output_dir="$DATADIR/$file_id"

# defaults
force=n
quiet=n
ic_file="$output_dir/$file_id.init"

# check for file_id
if [ -z $file_id ]; then
    echo "Error: file_id not set"
    echo "$usage"
    exit 1
fi

# shift the current values stored in the positional parameters
shift

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

# check if output_dir exists
if [ -d $output_dir ]; then
    if [ $force == "y" ]; then
	echo "Overwriting contents of $output_dir"
        rm -rf $output_dir/*
        mkdir $output_dir/images
    else
	echo "Error: $output_dir exists"
	echo "Use -f to override"
	exit 1
    fi
else
    mkdir $output_dir
    mkdir $output_dir/images
fi

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
sim_base="$sim_base --graph-dir $output_dir/images $sim_options"
sim_base="$sim_base $options"
sim_command="$sim_base --write-file $output_dir/$file_id"

# execute 
$sim_command

echo 

# averaging runs
echo "Averaging runs..."

avg_command="$CODEDIR/graph/bin/average_runs -d $dt -o $output_dir/$file_id"
avg_command="$avg_command $output_dir/$file_id*.sim.dat"

# execute average
$avg_command

# check if ic-file was generated
if [ -s $ic_file ]; then
    echo -ne # ic_file exists
else
    echo "Error: ic_file $ic_file does not exist"
    exit 1
fi

echo "Making plots"
plot_command="$CODEDIR/do_all/simplot.bash $file_id"
$plot_command

rm "$file_id.sim.tmp.gp"

# ghostview
#    gv $output_dir/$file_id.pairs.ps &
#    gv $output_dir/$file_id.singlets.ps &
