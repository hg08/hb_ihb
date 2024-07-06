load '~/bin/vdos_colorstyle.gnuplot'
set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output "Fig4.eps"
set encoding iso_8859_1

set size square 1.10,1.4

set multiplot
xmin=0.95
xmax=6.05

set origin 0, 0 
set xlabel "d (\305)" 
set xrange [xmin : xmax]
set xtics 1 
set mxtics 2
set ytics 0.2 
set mytics 2
set format y "%3.1f"

#====
#plot
#====
set size 1.1, 1.4
set border 1+2+4+8
set ylabel "k_{LC}, k_{IHB}  (ps^{-1})"
#set label 1 "b" right at graph 0.05, graph 0.95
set yrange [0.0: 0.8] 
set arrow 3 from 1.0,0.142 to 6,0.142 nohead ls 4 
set label 6 "k_{bulk}=0.14 ps^{-1}" left at graph 0.02, graph 0.10
set label 7 "0.14" at 0.25,0.14 textcolor rgb '#000000'
plot [xmin:xmax] \
     '../3_analyze/output/1_case1_kkprime.dat' u 1:2 with linespoints ls 20022 title "scenario 1  (LC)",\
     '../3_analyze/output/2_case2_kkprime.dat' u 1:2 with linespoints ls 10022 title "scenario 2 (IHB)"
     #'../3_analyze/output/1_case1_kkprime.dat' u 1:2:5 w yerrorbars ls 10022 notitle,\
     #'../3_analyze/output/2_case2_kkprime.dat' u 1:2:5 w yerrorbars ls 10022 notitle,\

unset multiplot
set term wxt
