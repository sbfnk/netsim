#!/bin/sh

for I in images/frame* 
do
  neato -Tps $I -o $I.ps
done


