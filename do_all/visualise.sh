#!/bin/sh

usage="Usage: visualise.sh file_id run_no"

if [ -z $2 ]; then
  echo "$usage"
  exit 1
fi

images_dir=$DATADIR/$1/images
run_no=$2
run_dir=`find $DATADIR/$1 -type d -name run\* | grep -E "run0*$run_no"`/images

if [ -z "$run_dir" ]; then
  echo "ERROR: run $run_no does not have directory in $images_dir"
  exit 1
fi

frames=`find $run_dir -type f -name frame\* | grep -E "frame[0-9]+$"`
if [ -z "$frames" ]; then
  png_frames=`find $run_dir -type f -name frame\* | grep -E "frame[0-9]+\.png$"`
  if [ -n "$png_frames" ]; then
    qiv -W 400 $run_dir/frame???.png &
    exit 0
  else 
    echo "ERROR: no frames found in $run_dir"
    exit 1
  fi
fi

for I in $frames
do
  echo $I
  neato -Tgif $I -o $I.gif
done

gthumb $run_dir &

exit 0

