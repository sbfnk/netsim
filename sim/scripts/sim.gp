set terminal postscript eps enhanced color dashed defaultplex "Helvetica" 24

set bmargin 1
set tmargin 1
set rmargin 1
set lmargin 1

set xlabel "time"
set ylabel "" 

set format y "%3.1e"

set style line 1 lt 1 lc rgb "blue"  lw 2
set style line 2 lt 1 lc rgb "red"   lw 2
set style line 3 lt 1 lc rgb "green" lw 2
set style line 4 lt 2 lc rgb "blue"  lw 3
set style line 5 lt 2 lc rgb "red"   lw 3
set style line 6 lt 2 lc rgb "green" lw 3

#####################################

set output 'sim.ps'
set title "simulation"
plot 'FILE_ID.sim.dat' u 1:2 ls 1 title 'S-' w l \
   , 'FILE_ID.sim.dat' u 1:3 ls 2 title 'I-' w l \
   , 'FILE_ID.sim.dat' u 1:4 ls 3 title 'R-' w l \
   , 'FILE_ID.sim.dat' u 1:5 ls 4 title 'S+' w l \
   , 'FILE_ID.sim.dat' u 1:6 ls 5 title 'I+' w l \
   , 'FILE_ID.sim.dat' u 1:7 ls 6 title 'R+' w l

#####################################

set output 'sim_sir.ps'
set title "simulation"
plot 'FILE_ID.mf.dat' u 1:($2+$5) ls 1 title 'S' w l \
   , 'FILE_ID.mf.dat' u 1:($3+$6) ls 2 title 'I' w l \
   , 'FILE_ID.mf.dat' u 1:($4+$7) ls 3 title 'R' w l


#####################################

#set output 'pairs1.ps'
#set title "Simulation"
#plot 'FILE_ID.sim.dat' u 1:8 title 'SSd' w l \
#   , 'FILE_ID.sim.dat' u 1:9 title 'SId' w l \
#   , 'FILE_ID.sim.dat' u 1:10 title 'SRd' w l \
#   , 'FILE_ID.sim.dat' u 1:14 title 'IId' w l \
#   , 'FILE_ID.sim.dat' u 1:15 title 'IRd' w l \
#   , 'FILE_ID.sim.dat' u 1:19 title 'RRd' w l
#
######################################
#
#set output 'pairs2.ps'
#set title "Simulation"
#plot 'FILE_ID.sim.dat' u 1:29 title 'SSi' w l \
#   , 'FILE_ID.sim.dat' u 1:30 title 'SIi' w l \
#   , 'FILE_ID.sim.dat' u 1:31 title 'SRi' w l \
#   , 'FILE_ID.sim.dat' u 1:35 title 'IIi' w l \
#   , 'FILE_ID.sim.dat' u 1:36 title 'IRi' w l \
#   , 'FILE_ID.sim.dat' u 1:40 title 'RRi' w l
#
######################################
#
#set output 'pairs3.ps'
#set title "Simulation"
#plot 'FILE_ID.sim.dat' u 1:23 title 'ssd' w l \
#   , 'FILE_ID.sim.dat' u 1:24 title 'sid' w l \
#   , 'FILE_ID.sim.dat' u 1:25 title 'srd' w l \
#   , 'FILE_ID.sim.dat' u 1:26 title 'iid' w l \
#   , 'FILE_ID.sim.dat' u 1:27 title 'ird' w l \
#   , 'FILE_ID.sim.dat' u 1:28 title 'rrd' w l
#
######################################
#
#set output 'pairs4.ps'
#set title "Simulation"
#plot 'FILE_ID.sim.dat' u 1:44 title 'ssi' w l \
#   , 'FILE_ID.sim.dat' u 1:45 title 'sii' w l \
#   , 'FILE_ID.sim.dat' u 1:46 title 'sri' w l \
#   , 'FILE_ID.sim.dat' u 1:47 title 'iii' w l \
#   , 'FILE_ID.sim.dat' u 1:48 title 'iri' w l \
#   , 'FILE_ID.sim.dat' u 1:49 title 'rri' w l
#
######################################
#
#set output 'pairs5.ps'
#set title "Simulation"
#plot 'FILE_ID.sim.dat' u 1:11 title 'Ssd' w l \
#   , 'FILE_ID.sim.dat' u 1:12 title 'Sid' w l \
#   , 'FILE_ID.sim.dat' u 1:13 title 'Srd' w l \
#   , 'FILE_ID.sim.dat' u 1:17 title 'Iid' w l \
#   , 'FILE_ID.sim.dat' u 1:18 title 'Ird' w l \
#   , 'FILE_ID.sim.dat' u 1:22 title 'Rrd' w l
#
######################################
#
#set output 'pairs6.ps'
#set title "Simulation"
#plot 'FILE_ID.sim.dat' u 1:32 title 'Ssi' w l \
#   , 'FILE_ID.sim.dat' u 1:33 title 'Sii' w l \
#   , 'FILE_ID.sim.dat' u 1:34 title 'Sri' w l \
#   , 'FILE_ID.sim.dat' u 1:38 title 'Iii' w l \
#   , 'FILE_ID.sim.dat' u 1:39 title 'Iri' w l \
#   , 'FILE_ID.sim.dat' u 1:43 title 'Rri' w l
#
######################################
#
#set output 'pairs7.ps'
#set title "Simulation"
#plot 'FILE_ID.sim.dat' u 1:16 title 'sId' w l \
#   , 'FILE_ID.sim.dat' u 1:20 title 'sRd' w l \
#   , 'FILE_ID.sim.dat' u 1:21 title 'iRd' w l \
#
######################################
#
#set output 'pairs8.ps'
#set title "Simulation"
#plot 'FILE_ID.sim.dat' u 1:37 title 'sIi' w l \
#   , 'FILE_ID.sim.dat' u 1:41 title 'sRi' w l \
#   , 'FILE_ID.sim.dat' u 1:42 title 'iRi' w l \
#
######################################
#
#### CHANGE AS NEEDED !!!!!!!!!!!!!
##N=1e5
##Qd=3
##Qi=3
#### CHANGE AS NEEDED !!!!!!!!!!!!!
#
#set output 'Cxx1.ps'
#set title "Correlation"
#set xrange [0:20]
#set yrange [0:1]
#set key spacing 2

# try1:
#plot 'FILE_ID.sim.dat' u 1:(N*$9/($2*$3)/Qd) title 'CSId' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$30/($2*$3)/Qi) title 'CSIi' w l \
#   , 'FILE_ID.sim.dat' u 1:($3/N) title 'I' w l \

# try2:
#plot 'FILE_ID.sim.dat' u 1:(N*$20/($2*$5)/Qd) title 'CSsd' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$41/($2*$5)/Qi) title 'CSsi' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$29/($2*$2)/Qi) title 'CSSi' w l \

# try3:
#plot 'FILE_ID.sim.dat' u 1:(N*$20/($2*$5)/Qd) title 'CSsd' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$41/($2*$5)/Qi) title 'CSsi' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$29/($2*$2)/Qi) title 'CSSi' w l \
#   , 'FILE_ID.sim.dat' u 1:($5/N) title 's' w l \

# try4:
#plot 'FILE_ID.sim.dat' u 1:(N*$26/($5*$3)/Qd) title 'CsId' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$9/($2*$3)/Qd) title 'CSId' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$41/($2*$5)/Qi) title 'CSsi' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$29/($2*$2)/Qi) title 'CSSi' w l \
#   , 'FILE_ID.sim.dat' u 1:($3/N) title 'I' w l \

# try5
#plot 'FILE_ID.sim.dat' u 1:(N*$15/($5*$6)/Qd) title 'Csid' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$36/($5*$6)/Qi) title 'Csii' w l \
#   , 'FILE_ID.sim.dat' u 1:($5/N) title 's' w l \
#   , 'FILE_ID.sim.dat' u 1:($6/N) title 'i' w l \

# try6
#plot 'FILE_ID.sim.dat' u 1:(N*$15/($5*$6)/Qd) title 'Csid' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$36/($5*$6)/Qi) title 'Csii' w l \
#   , 'FILE_ID.sim.dat' u 1:($3/N) title 'I' w l \
#   , 'FILE_ID.sim.dat' u 1:($5/N) title 's' w l \
#   , 'FILE_ID.sim.dat' u 1:($6/N) title 'i' w l \

# try7
#plot 'FILE_ID.sim.dat' u 1:(N*$15/($5*$6)/Qd) title 'Csid' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$26/($5*$3)/Qd) title 'CsId' w l \
#   , 'FILE_ID.sim.dat' u 1:(($3+$6)/N) title 'I+i' w l \

# try10
#plot 'FILE_ID.sim.dat' u 1:(N*$9/($2*$3)/Qd) title 'CSId' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$21/($2*$6)/Qd) title 'CSid' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$26/($5*$3)/Qd) title 'CsId' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$15/($5*$6)/Qd) title 'Csid' w l \
#   , 'FILE_ID.sim.dat' u 1:(($3+$6)/N) title 'I+i' w l \

#set output 'Cxx2.ps'
#set title "Correlation"
#set xrange [0:20]
#set yrange [0:2.5]
#set key spacing 2

#plot 'FILE_ID.sim.dat' u 1:(N*$30/($2*$3)/Qi) title 'CSIi' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$42/($2*$6)/Qi) title 'CSii' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$47/($5*$3)/Qi) title 'CsIi' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$36/($5*$6)/Qi) title 'Csii' w l \
#   , 'FILE_ID.sim.dat' u 1:(($3+$6)/N) title 'I+i' w l \

#set output 'Cxx3.ps'
#set title "Correlation"
#set xrange [0:20]
#set yrange [0:2.5]
#set key spacing 2

#plot 'FILE_ID.sim.dat' u 1:(N*$8/($2*$2)/Qd) title 'CSSd' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$20/($2*$5)/Qd) title 'CSsd' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$14/($5*$5)/Qd) title 'Cssd' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$29/($2*$2)/Qi) title 'CSSi' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$41/($2*$5)/Qi) title 'CSsi' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$35/($5*$5)/Qi) title 'Cssi' w l \

#set output 'Cxx4.ps'
#set title "Correlation"
#set xrange [0:20]
#set yrange [0:2.5]
#set key spacing 2

#plot 'FILE_ID.sim.dat' u 1:(N*$11/($3*$3)/Qd) title 'CIId' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$23/($3*$6)/Qd) title 'CIid' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$17/($6*$6)/Qd) title 'Ciid' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$32/($3*$3)/Qi) title 'CIIi' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$44/($3*$6)/Qi) title 'CIii' w l \
#   , 'FILE_ID.sim.dat' u 1:(N*$38/($6*$6)/Qi) title 'Ciii' w l \






quit



