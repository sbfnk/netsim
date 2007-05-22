#!/bin/sh

for I in images/frame* 
do
  neato -Tps2 $I -o $I.ps
done


