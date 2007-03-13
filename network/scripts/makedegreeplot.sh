#!/bin/bash

# set script path
gp_script=degree.gp
latex_script=fix.tex
join_script=join.tex

# check command-line args
if [[ $# != 1 ]]; then
    echo "Usage: ./makedegreeplot.bash file_id"
    echo "       file_id is missing."
    exit 1
fi

# create gnuplot script for ode
sed -e s:FILE_ID:$1:g $gp_script > $1.tmp.gp

# make figs
gnuplot $1.tmp.gp > /dev/null
rm -f $1.tmp.gp

# ps -> eps
ps2epsi degree_dist.ps degree_dist.eps
ps2epsi degree_dist_log.ps degree_dist_log.eps

# create latex script
sed -e s:FILE_ID:degree_dist:g $latex_script > $1.tmp.tex

# fix fonts
latex $1.tmp.tex > /dev/null 
dvips -q -o $1.degree.ps $1.tmp.dvi > /dev/null
rm -f $1.tmp.log $1.tmp.aux $1.tmp.dvi $1.tmp.tex

# create latex script for log-log plot
sed -e s:FILE_ID:degree_dist_log:g $latex_script > $1.tmp.tex

# fix fonts
latex $1.tmp.tex > /dev/null 
dvips -q -o $1.degree_log.ps $1.tmp.dvi > /dev/null
rm -f $1.tmp.log $1.tmp.aux $1.tmp.dvi $1.tmp.tex

# fix lines
sed -e "s/PL \[4 dl1 2 dl2\]/PL \[12 dl1 6 dl2\]/g" $1.degree.ps > tmp.ps
mv -f tmp.ps $1.degree.ps
ps2epsi $1.degree.ps $1.degree.eps

sed -e "s/PL \[4 dl1 2 dl2\]/PL \[12 dl1 6 dl2\]/g" $1.degree_log.ps > tmp.ps
mv -f tmp.ps $1.degree_log.ps
ps2epsi $1.degree_log.ps $1.degree_log.eps

# join
sed -e s:FILE_ID:$1.degree:g $join_script > $1.tmp.tex
latex $1.tmp.tex > /dev/null
dvips -q -o $1.degrees.ps $1.tmp.dvi > /dev/null
rm -f $1.tmp.log $1.tmp.aux $1.tmp.dvi $1.tmp.tex

ps2epsi $1.degrees.ps $1.degrees.eps


# clean
rm -f degree_dist.ps degree_dist.eps  $1.degree.ps $1.degree_log.ps
rm -f  $1.degree.eps $1.degree_log.eps $1.degrees.ps
rm -f degree_dist_log.ps degree_dist_log.eps

# gv
gv $1.degrees.eps &

exit 0
