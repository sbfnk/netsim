#!/bin/sh
# script for animation of gif images
# usage: mkanim [delay in 1/100 seconds]
# default delay is specified below

ANIM_CMD=gifsicle
DELAY=10 # 0.1 seconds

program=$(which $ANIM_CMD)
    
if [ -z "$program" ] ; then
	echo "Program $ANIM_CMD could not be found."
	exit -1
fi

if [ -n "$1" ] ; then
  ANIM_OPT="--delay $1"
else 
  ANIM_OPT="--delay $DELAY"
fi

echo "creating animated gif..."
$ANIM_CMD $ANIM_OPT --no-loopcount --optimize $(ls images/frame*.gif) > animation.gif

echo "...saved as animation.gif"

exit 0
