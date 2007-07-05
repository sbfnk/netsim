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
usage="Usage: do_all.sh file_id [-m model_file] [-o ode_file]"
usage="$usage [-s sim_file] [-d timestep] [-q] [model_options]"

# set file_id
file_id=$1

# set output dir
output_dir="$DATADIR/$file_id"

# defaults
force=n
quiet=n
ode_solver='solve_ode.x'
ic_file="$output_dir/$file_id.init"

# check for file_id
if [ -z $file_id ]; then
    echo "Error: file_id not set"
    echo "$usage"
    exit 1
fi

# saving command line 
comm_line="$0 $@"

# shift the current values stored in the positional parameters
shift

# parse argument list
while [ $# -ge 1 ]; 
  do
  
  case $1 
      in
      -m) shift; model=$1;;
      -o) shift; ode=$1;;
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

# writing command line
echo $comm_line > $output_dir/$file_id.comm_line

# setting sim_options
sim_options=""

if [ $model ]; then
    sim_options="$sim_options -m $model"
fi
if [ $sim ]; then
    sim_options="$sim_options -p $sim"
fi

# saving parameters files
cp -f $model $output_dir/$file_id.model.prm
cp -f $sim $output_dir/$file_id.sim.prm

# set sim_command
sim_base="$CODEDIR/graph/bin/simulate"
sim_base="$sim_base --graph-dir $output_dir/images $sim_options"
sim_base="$sim_base $options"
sim_command="$sim_base --write-file $output_dir/$file_id

# execute 
$sim_command

# averaging runs
echo
echo "Averaging runs..."

avg_command="$CODEDIR/graph/bin/average_runs -d $dt -o $output_dir/$file_id"
avg_command="$avg_command $output_dir/$file_id*.sim.dat"

# execute average
$avg_command

rm $output_dir/$file_id.sim.tmp.gp

# check if ic-file was generated
if [ -s $ic_file ]; then
    echo -ne # ic_file exists
else
    echo "Error: ic_file $ic_file does not exist"
    exit 1
fi

echo "Running ode..."

# set ode/model prm files
if [ $model ]; then
    ode_options="$ode_options -m $model"
fi
if [ $ode ]; then
    ode_options="$ode_options -p $ode"
fi

# copy ode params file
cp -f $ode $output_dir/$file_id.ode.prm

# set ode command
ode_command="$CODEDIR/ode_solver/$ode_solver --file-id $output_dir/$file_id"
ode_command="$ode_command --ic-file $output_dir/$file_id.init"
ode_command="$ode_command $ode_options --di-model --mf-model $options"

# execute ode solver
$ode_command

# plot
if [ $quiet == "n" ]; then
    echo "Making plots"
    plot_command="$CODEDIR/do_all/makeplot.bash $file_id"
    $plot_command
    
# ghostview
#    gv $output_dir/$file_id.pairs.ps &
#    gv $output_dir/$file_id.singlets.ps &
fi
