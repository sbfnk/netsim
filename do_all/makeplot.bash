#!/bin/bash

# set script path
ode_dir=$CODEDIR/ode_solver
graph_dir=$CODEDIR/graph

gp_ode_script=$ode_dir/mfpa.gp
gp_sim_script=$graph_dir/scripts/sim.gp

# check command-line args
if [[ $# != 1 ]]; then
    echo "Usage: ./makeplot.bash file_id"
    echo "       file_id is missing."
    exit 1
fi

# output directory
output_dir=$DATADIR/$1

# create gnuplot script for ode
cat $output_dir/$1.gp > $output_dir/$1.ode.tmp.gp
sed -e s:FILE_ID:$output_dir/$1:g $gp_ode_script >> $output_dir/$1.ode.tmp.gp

# make figs
gnuplot $output_dir/$1.ode.tmp.gp > /dev/null
rm -f $output_dir/$1.ode.tmp.gp

# ps -> eps
ps2epsi mf.ps mf.eps
ps2epsi pa-di.ps pa-di.eps
ps2epsi pa-dib.ps pa-dib.eps
#ps2epsi pa.ps pa.eps
#ps2epsi pa1.ps pa1.eps
#ps2epsi pa2.ps pa2.eps
#ps2epsi pa3.ps pa3.eps
#ps2epsi pa4.ps pa4.eps
#ps2epsi pa5.ps pa5.eps
#ps2epsi pa6.ps pa6.eps
#ps2epsi pa7.ps pa7.eps
#ps2epsi pa8.ps pa8.eps
#ps2epsi Cxx1.ps Cxx1.eps
#ps2epsi Cxx2.ps Cxx2.eps
#ps2epsi Cxx3.ps Cxx3.eps
#ps2epsi Cxx4.ps Cxx4.eps

# fix fonts
cp -f $ode_dir/fix1.tex .
latex fix1.tex > /dev/null 
dvips -q -o $output_dir/$1.mfpa.ps fix1.dvi > /dev/null
rm -f fix1.tex fix1.log fix1.aux fix1.dvi

# fix lines
sed -e "s/PL \[4 dl1 2 dl2\]/PL \[12 dl1 6 dl2\]/g" $output_dir/$1.mfpa.ps > tmp.ps
mv -f tmp.ps $output_dir/$1.mfpa.ps

# clean
rm -f mf.ps mf.eps
rm -f pa-di.ps pa-di.eps
rm -f pa-dib.ps pa-dib.eps

# create gnuplot script for simulaion
cat $output_dir/$1.gp > $output_dir/$1.sim.tmp.gp
sed -e s:FILE_ID:$output_dir/$1:g $gp_sim_script >> $output_dir/$1.sim.tmp.gp

# make figs
gnuplot $output_dir/$1.sim.tmp.gp > /dev/null
rm -f $output_dir/$1.sim.tmp.gp

# ps -> eps
ps2epsi sim.ps sim.eps
ps2epsi mf.ps mf.eps
ps2epsi pa-di.ps pa-di.eps
ps2epsi pa-dib.ps pa-dib.eps
#ps2epsi sim_sir.ps sim_sir.eps
#ps2epsi pairs1.ps pairs1.eps
#ps2epsi pairs2.ps pairs2.eps
#ps2epsi pairs3.ps pairs3.eps
#ps2epsi pairs4.ps pairs4.eps
#ps2epsi pairs5.ps pairs5.eps
#ps2epsi pairs6.ps pairs6.eps
#ps2epsi pairs7.ps pairs7.eps
#ps2epsi pairs8.ps pairs8.eps

# create ps files
latex singlets.tex > /dev/null 
dvips -q -o $output_dir/$1.singlets.ps singlets.dvi > /dev/null
#latex pairs.tex  > /dev/null 
#dvips -q -o $output_dir/$1.pairs.ps pairs.dvi > /dev/null

# fix lines
sed -e "s/PL \[4 dl1 2 dl2\]/PL \[12 dl1 6 dl2\]/g" $output_dir/$1.singlets.ps > tmp.ps
mv -f tmp.ps $output_dir/$1.singlets.ps

mv sim.eps $output_dir/$1.sim.eps
mv mf.eps $output_dir/$1.mf.eps
mv pa-di.eps $output_dir/$1.pa-di.eps
mv pa-dib.eps $output_dir/$1.pa-dib.eps
#mv sim_sir.eps $output_dir/$1.sim_sir.eps
#mv mf_sir.eps $output_dir/$1.mf_sir.eps

# clean
rm -f sim.ps sim.eps
rm -f mf.ps mf.eps
rm -f pa-di.ps pa-di.eps
rm -f pa-dib.ps pa-dib.eps
rm -f singlets.log singlets.aux singlets.dvi 
rm -f pairs.log pairs.aux pairs.dvi 
rm -f mf.ps mf_sir.ps
rm -f pa-di.ps pa-di.eps
rm -f pa-dib.ps pa-dib.eps
rm -f pa1.ps pa1.eps 
rm -f pa2.ps pa2.eps 
rm -f pa3.ps pa3.eps 
rm -f pa4.ps pa4.eps
rm -f pa5.ps pa5.eps 
rm -f pa6.ps pa6.eps 
rm -f pa7.ps pa7.eps 
rm -f pa8.ps pa8.eps
rm -f Cxx1.ps Cxx1.eps
rm -f Cxx2.ps Cxx2.eps
rm -f Cxx3.ps Cxx3.eps
rm -f Cxx4.ps Cxx4.eps
rm -f sim.ps sim_sir.ps
rm -f pairs1.ps pairs1.eps
rm -f pairs2.ps pairs2.eps
rm -f pairs3.ps pairs3.eps
rm -f pairs4.ps pairs4.eps
rm -f pairs5.ps pairs5.eps
rm -f pairs6.ps pairs6.eps
rm -f pairs7.ps pairs7.eps
rm -f pairs8.ps pairs8.eps

# extracting single figures
#for fig in mf.eps pa.eps pa5.eps pa6.eps Cxx.eps
#do
#    echo "... $fig"
#    sed -e s/FIGURE/$fig/g $fixf > tmp.tex
#    latex tmp.tex > /dev/null  
#    dvips -q -o tmp.ps tmp.dvi > /dev/null
#    ps2epsi tmp.ps $1.$fig
#    rm -f tmp.dvi tmp.log tmp.aux tmp.ps tmp.tex $fig
#done	

exit 0
