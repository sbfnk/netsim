set terminal postscript eps enhanced defaultplex \
    color dashed dashlength 3.0 linewidth 1. \
    "Helvetica" 24 

set bmargin 1
set tmargin 1
set rmargin 1
set lmargin 1

set xlabel "PARAM";
set ylabel "";

set format y "%3.1e"

set style line 1 lt 1 lc rgb "blue"  lw 2
set style line 2 lt 1 lc rgb "red"   lw 2
set style line 3 lt 1 lc rgb "green" lw 2
set style line 4 lt 2 lc rgb "blue"  lw 3
set style line 5 lt 2 lc rgb "red"   lw 3
set style line 6 lt 2 lc rgb "green" lw 3

#####################################

set output 'mf.ps'
set title "Mean Field"
plot 'FILE_ID.mf.dat' u 1:2 ls 1 title 'S-' w l \
   , 'FILE_ID.mf.dat' u 1:3 ls 2 title 'I-' w l \
   , 'FILE_ID.mf.dat' u 1:4 ls 3 title 'R-' w l \
   , 'FILE_ID.mf.dat' u 1:5 ls 4 title 'S+' w l \
   , 'FILE_ID.mf.dat' u 1:6 ls 5 title 'I+' w l \
   , 'FILE_ID.mf.dat' u 1:7 ls 6 title 'R+' w l

#####################################

#set output 'pa.ps'
#set title "Pair Approx"
#plot 'FILE_ID.pa.dat' u 1:2 ls 1 title 'S-' w l \
#   , 'FILE_ID.pa.dat' u 1:3 ls 2 title 'I-' w l \
#   , 'FILE_ID.pa.dat' u 1:4 ls 3 title 'R-' w l \
#   , 'FILE_ID.pa.dat' u 1:5 ls 4 title 'S+' w l \
#   , 'FILE_ID.pa.dat' u 1:6 ls 5 title 'I+' w l \
#   , 'FILE_ID.pa.dat' u 1:7 ls 6 title 'R+' w l

#####################################

set output 'sim.ps'
set title "Simluation"
plot 'FILE_ID.sim.dat' u 1:2 ls 1 title 'S-' w l \
   , 'FILE_ID.sim.dat' u 1:3 ls 2 title 'I-' w l \
   , 'FILE_ID.sim.dat' u 1:4 ls 3 title 'R-' w l \
   , 'FILE_ID.sim.dat' u 1:5 ls 4 title 'S+' w l \
   , 'FILE_ID.sim.dat' u 1:6 ls 5 title 'I+' w l \
   , 'FILE_ID.sim.dat' u 1:7 ls 6 title 'R+' w l

quit



