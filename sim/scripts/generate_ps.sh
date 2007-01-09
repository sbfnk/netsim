#!/bin/sh

neato -Gstart=rand -Tps images/start -o images/a_start.ps
for I in images/frame* 
do
  neato -Gstart=rand -Tps $I -o $I.ps
done
neato -Gstart=rand -Tps images/end -o images/z_end.ps


