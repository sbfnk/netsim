#!/bin/sh
# script for animation of gif images
# usage: mkanim [delay in 1/100 seconds]
# default delay is specified below

#ANIM_CMD=gifsicle
ANIM_CMD=animate
#DISP_CMD=firefox
DELAY=10 # 0.1 seconds

program=$(which $ANIM_CMD)
    
usage="Usage: anmiate.sh file_id run_no [anim_delay]"

if [ -z $2 ]; then
  echo "$usage"
  exit 1
fi

images_dir=$DATADIR/$1/images
run_no=$2
run_dir=`find $DATADIR/$1/images -type d -name run\* | grep -E "run0*$run_no"`

#rm $images_dir/animation.gif

if [ -n "$3" ] ; then
  DELAY=$3
fi

if [ -z "$run_dir" ]; then
  echo "ERROR: run $run_no does not have directory in $images_dir"
  exit 1
fi

gif_frames=`find $run_dir -type f -name frame\* | grep -E "frame[0-9]+\.gif$"`
if [ -z "$gif_frames" ]; then
  png_frames=`find $run_dir -type f -name frame\* | grep -E "frame[0-9]+\.png$"`
  if [ -n "$png_frames" ]; then
    qiv -s -d $(echo "scale=1;$DELAY/100" | bc) -W 400 $run_dir/frame???.png &
    exit 0
  else 
    echo "ERROR: no frames found in $run_dir"
    exit 1
  fi
fi

$ANIM_CMD --delay $DELAY $run_dir/frame???.gif &

#echo "creating animated gif..."
#$ANIM_CMD $ANIM_OPT --no-loopcount --optimize $(ls $images_dir/frame*.gif) > $images_dir/animation.gif

#echo "...saved as $images_dir/animation.gif"

#animate=$(which $DISP_CMD)

#if [ -z "$animate" ]; then
#  echo "$DISP_CMD not found, will not play animation"
#else
#  $animate $images_dir/animation.gif
#fi

exit 0
