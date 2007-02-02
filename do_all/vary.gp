set terminal postscript eps enhanced color dashed defaultplex "Helvetica" 24

set bmargin 1
set tmargin 1
set rmargin 1
set lmargin 1

set xlabel "PARAM" 0.000000,0.000000
set ylabel "" 0.000000,0.000000

set format y "%3.1e"

#####################################

set output 'mf.ps'
set title "Mean Field"
plot 'FILE_ID.mf.dat' u 1:2 title 'S-' w l \
   , 'FILE_ID.mf.dat' u 1:3 title 'I-' w l \
   , 'FILE_ID.mf.dat' u 1:4 title 'R-' w l \
   , 'FILE_ID.mf.dat' u 1:5 title 'S+' w l \
   , 'FILE_ID.mf.dat' u 1:6 title 'I+' w l \
   , 'FILE_ID.mf.dat' u 1:7 title 'R+' w l

#####################################

set output 'pa.ps'
set title "Pair Approx"
plot 'FILE_ID.pa.dat' u 1:2 title 'S-' w l \
   , 'FILE_ID.pa.dat' u 1:3 title 'I-' w l \
   , 'FILE_ID.pa.dat' u 1:4 title 'R-' w l \
   , 'FILE_ID.pa.dat' u 1:5 title 'S+' w l \
   , 'FILE_ID.pa.dat' u 1:6 title 'I+' w l \
   , 'FILE_ID.pa.dat' u 1:7 title 'R+' w l

#####################################

set output 'sim.ps'
set title "Simluation"
plot 'FILE_ID.sim.dat' u 1:2 title 'S-' w l \
   , 'FILE_ID.sim.dat' u 1:3 title 'I-' w l \
   , 'FILE_ID.sim.dat' u 1:4 title 'R-' w l \
   , 'FILE_ID.sim.dat' u 1:5 title 'S+' w l \
   , 'FILE_ID.sim.dat' u 1:6 title 'I+' w l \
   , 'FILE_ID.sim.dat' u 1:7 title 'R+' w l

quit



