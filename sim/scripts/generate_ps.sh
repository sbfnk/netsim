#!/bin/sh

neato -Tps images/start -o images/a_start.ps
for I in images/frame* 
do
  neato -Tps $I -o $I.ps
done
neato -Tps images/end -o images/z_end.ps


