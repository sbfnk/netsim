#!/bin/bash

# set script path
graph_dir=$CODEDIR/graph

degree_script=$graph_dir/scripts/makedegreeplot.sh

# check command-line args
if [[ $# != 1 ]]; then
    echo "Usage: ./makedegree.sh file_id"
    echo "       file_id is missing."
    exit 1
fi

# output directory
output_dir=$DATADIR/$1

# make degree plot
if [[ -f $output_dir/$1.degree ]]; then
  $degree_script $output_dir/$1
else
  echo "ERROR: no degree file $output_dir/$1.degree"
fi

exit 0
