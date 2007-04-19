#!/bin/sh

rm images/*.gif
for I in images/frame* 
do
  echo $I
  neato -Tgif $I -o $I.gif
done


