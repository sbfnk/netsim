#!/bin/sh

rm images/*.gif
neato -Tgif images/start -o images/a_start.gif
for I in images/frame* 
do
  echo $I
  neato -Tgif $I -o $I.gif
done
neato -Tgif images/end -o images/z_end.gif
cd -


