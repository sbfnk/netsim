set terminal postscript eps enhanced color dashed defaultplex "Helvetica" 24

set bmargin 1
set tmargin 1
set rmargin 1
set lmargin 1

set xlabel "time"
set ylabel "" 

set format y "%3.1e"

set key right top samplen 2 spacing 1.5
set xrange [0:Tmax]
set yrange [0:N]

set style line 1 lt 1 lc rgb "blue"  lw 2
set style line 2 lt 2 lc rgb "blue"  lw 3
set style line 3 lt 1 lc rgb "red"   lw 2
set style line 4 lt 1 lc rgb "green" lw 2

#####################################

set output 'sim.ps'
set title "simulation"
plot 'FILE_ID.sim.dat' u 1:2 ls 1 title 'S-' w l \
   , 'FILE_ID.sim.dat' u 1:3 ls 2 title 'S+' w l \
   , 'FILE_ID.sim.dat' u 1:4 ls 3 title 'I' w l \
   , 'FILE_ID.sim.dat' u 1:5 ls 4 title 'R' w l \

#####################################

#set output 'mf.ps'
#set title "Mean Field"
#plot 'FILE_ID.mf.dat' u 1:2 title 'S-' w l ls 1 \
#   , 'FILE_ID.mf.dat' u 1:3 title 'I-' w l ls 2 \
#   , 'FILE_ID.mf.dat' u 1:4 title 'R-' w l ls 3 \
#   , 'FILE_ID.mf.dat' u 1:5 title 'S+' w l ls 4 \
#   , 'FILE_ID.mf.dat' u 1:6 title 'I+' w l ls 5 \
#   , 'FILE_ID.mf.dat' u 1:7 title 'R+' w l ls 6

#####################################

#set output 'pa-di.ps'
#set title "Pair Approx di"
#plot 'FILE_ID.di.dat' u 1:2 title 'S-' w l ls 1 \
#   , 'FILE_ID.di.dat' u 1:3 title 'I-' w l ls 2 \
#   , 'FILE_ID.di.dat' u 1:4 title 'R-' w l ls 3 \
#   , 'FILE_ID.di.dat' u 1:5 title 'S+' w l ls 4 \
#   , 'FILE_ID.di.dat' u 1:6 title 'I+' w l ls 5 \
#   , 'FILE_ID.di.dat' u 1:7 title 'R+' w l ls 6

#####################################

#set output 'pa-dib.ps'
#set title "Pair Approx dib"
#plot 'FILE_ID.dib.dat' u 1:2 title 'S-' w l ls 1 \
#   , 'FILE_ID.dib.dat' u 1:3 title 'I-' w l ls 2 \
#   , 'FILE_ID.dib.dat' u 1:4 title 'R-' w l ls 3 \
#   , 'FILE_ID.dib.dat' u 1:5 title 'S+' w l ls 4 \
#   , 'FILE_ID.dib.dat' u 1:6 title 'I+' w l ls 5 \
#   , 'FILE_ID.dib.dat' u 1:7 title 'R+' w l ls 6

quit

