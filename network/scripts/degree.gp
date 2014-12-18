set terminal postscript eps enhanced color dashed defaultplex "Helvetica" 24

set bmargin 1
set tmargin 1
set rmargin 1
set lmargin 1

set xlabel "degree"
set ylabel "" 


set format y "%3.1e"

set key right top samplen 2.5 spacing 1.3
#set xrange [0:MaxDegree]
#set yrange [0:N]

set style line 1 lt 2 lc "black" lw 2
set style line 2 lt 3 lc "black" lw 2
set style line 3 lt 1 lc "black" lw 2

#####################################

set output 'degree_dist.ps'
set title "degdist"
plot 'FILE_ID.degree' u 1:2 ls 1 title 'dedge---' w lp \
   , 'FILE_ID.degree' u 1:3 ls 2 title 'iedge---' w lp \
   , 'FILE_ID.degree' u 1:4 ls 3 title 'total---' w lp 

#####################################

set logscale xy

set output 'degree_dist_log.ps'
set title "degdist"
plot 'FILE_ID.degree' u 1:2 ls 1 title 'dedge---' w lp \
   , 'FILE_ID.degree' u 1:3 ls 2 title 'iedge---' w lp \
   , 'FILE_ID.degree' u 1:4 ls 3 title 'total---' w lp 


quit



