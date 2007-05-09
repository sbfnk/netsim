#!/bin/bash

# set script path
gp_script=$CODEDIR/graph/scripts/degree.gp
latex_script=$CODEDIR/graph/scripts/fix.tex
join_script=$CODEDIR/graph/scripts/join.tex

# check command-line args
if [[ $# != 1 ]]; then
    echo "Usage: ./makedegreeplot.bash file_id"
    echo "       file_id is missing."
    exit 1
fi

# create gnuplot script for ode
sed -e s:FILE_ID:$1:g $gp_script > degree.tmp.gp

# make figs
gnuplot degree.tmp.gp > /dev/null
#rm -f degree.tmp.gp

# ps -> eps
ps2epsi degree_dist.ps degree_dist.eps
ps2epsi degree_dist_log.ps degree_dist_log.eps

# create latex script
sed -e s:FILE_ID:degree_dist:g $latex_script > degree.tmp.tex

# fix fonts
latex degree.tmp.tex > /dev/null 
dvips -q -o degree.ps degree.tmp.dvi > /dev/null
#rm -f degree.tmp.log degree.tmp.aux degree.tmp.dvi degree.tmp.tex

# create latex script for log-log plot
sed -e s:FILE_ID:degree_dist_log:g $latex_script > degree.tmp.tex

# fix fonts
latex degree.tmp.tex > /dev/null 
dvips -q -o degree_log.ps degree.tmp.dvi > /dev/null
#rm -f degree.tmp.log degree.tmp.aux degree.tmp.dvi degree.tmp.tex

# fix lines
sed -e "s/PL \[4 dl1 2 dl2\]/PL \[12 dl1 6 dl2\]/g" degree.ps > tmp.ps
mv -f tmp.ps degree.ps 
ps2epsi degree.ps degree_fixed.eps

sed -e "s/PL \[4 dl1 2 dl2\]/PL \[12 dl1 6 dl2\]/g" degree_log.ps > tmp.ps
mv -f tmp.ps degree_log.ps
ps2epsi degree_log.ps degree_fixed_log.eps

# join
sed -e s:FILE_ID:degree_fixed:g $join_script > degree.tmp.tex
latex degree.tmp.tex > /dev/null
dvips -q -o degrees.ps degree.tmp.dvi > /dev/null
#rm -f $1.tmp.log $1.tmp.aux $1.tmp.dvi $1.tmp.tex

mv degree_fixed.eps $1.degree.eps
mv degree_fixed_log.eps $1.degree_log.eps
ps2epsi degrees.ps $1.degrees.eps

# clean
rm -f degree_dist.eps degree_dist_log.eps
rm -f degree_dist.ps degree_dist_log.ps
rm -f degree.ps degree_log.ps degrees.ps
rm -f degree.tmp.*

# gv
gv $1.degrees.eps &

exit 0
