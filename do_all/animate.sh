#!/bin/sh
# script for animation of gif images
# usage: mkanim [delay in 1/100 seconds]
# default delay is specified below

ANIM_CMD=gifsicle
DISP_CMD=firefox
DELAY=10 # 0.1 seconds

program=$(which $ANIM_CMD)
    
if [ -z "$program" ] ; then
	echo "Program $ANIM_CMD could not be found."
	exit -1
fi

usage="Usage: anmiate.sh file_id [anim_delay]"

if [ -z $1 ]; then
  echo "Error: file_id not set"
  echo "$usage"
  exit 1
fi

images_dir=$DATADIR/$1/images

rm $images_dir/animation.gif

if [ -n "$2" ] ; then
  ANIM_OPT="--delay $2"
else 
  ANIM_OPT="--delay $DELAY"
fi

echo "creating animated gif..."
$ANIM_CMD $ANIM_OPT --no-loopcount --optimize $(ls $images_dir/frame*.gif) > $images_dir/animation.gif

echo "...saved as $images_dir/animation.gif"

animate=$(which $DISP_CMD)

if [ -z "$animate" ]; then
  echo "$DISP_CMD not found, will not play animation"
else
  $animate $images_dir/animation.gif
fi

exit 0
