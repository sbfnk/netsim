#!/bin/sh

usage="Usage: visualise.sh file_id"

if [ -z $1 ]; then
  echo "Error: file_id not set"
  echo "$usage"
  exit 1
fi

images_dir=$DATADIR/$1/images

sed -i 's/^.*rewired.*$//g' $images_dir/*

rm $images_dir/*.gif
for I in $images_dir/frame* 
do
  echo $I
  neato -Tgif $I -o $I.gif
done

gthumb $images_dir &
