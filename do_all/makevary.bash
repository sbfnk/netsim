#!/bin/bash

gp_script=$CODEDIR/do_all/vary.gp
fixf=$CODEDIR/ode/fix21.tex

# check command-line args
if [[ $# != 2 ]]; then
    echo "Usage: ./makevary.bash file_id param"
    echo "       one argument is missing."
    exit 1
fi

output_dir=$DATADIR/$1

# create gnuplot script for ode
cp $output_dir/$1.gp $output_dir/$1.tmp.gp
cat $gp_script >> $output_dir/$1.tmp.gp
sed -i s:FILE_ID:$output_dir/$1:g $output_dir/$1.tmp.gp
sed -i s:PARAM:$2:g $output_dir/$1.tmp.gp

# make figs
gnuplot $output_dir/$1.tmp.gp > /dev/null
#rm -f $output_dir/$1.tmp.gp

# ps -> eps
ps2epsi mf.ps mf.eps
ps2epsi pa.ps pa.eps
ps2epsi sim.ps sim.eps

# create ps files
latex vary.tex > /dev/null 
dvips -q -o $output_dir/$1.ps vary.dvi > /dev/null

# clean
rm -f vary.log vary.aux vary.dvi 
rm -f mf.ps pa.ps sim.ps
