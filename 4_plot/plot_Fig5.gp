load 'colorstyle.gnuplot'
set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output "Fig5.eps" 
set encoding iso_8859_1

set size square 2.10,1.4

set multiplot
# Plot c(t*) of d
xmin=0.95
xmax=6.05

#====
#plot
#====
set origin 0, 0 
set size 1.1, 1.4
set xlabel "d (\305)" 
set border 1+2+4+8
set xrange [xmin : xmax]
set yrange [0 : 1.0]
set xtics 1 
set ytics 0.2 
set format y "%2.1f"
set ylabel "C_2(t*)"
set key at graph 0.6, graph 0.9
set label 1 "(a)" left at graph 0.01, graph 0.95



plot \
    '../3_analyze/output/layers_c2_at_ref.dat' i 0 u 1:2 w l ls 10021 title "t* = 1 ps",\
    '../3_analyze/output/layers_c2_at_ref.dat' i 0 u 1:2:3 w yerrorbars ls 10021 notitle,\
    '../3_analyze/output/layers_c2_at_ref.dat' i 1 u 1:2 w l ls 10022 title "t* = 2 ps",\
    '../3_analyze/output/layers_c2_at_ref.dat' i 1 u 1:2:3 w yerrorbars ls 10022 notitle,\
    '../3_analyze/output/layers_c2_at_ref.dat' i 2 u 1:2 w l ls 10023 title "t* = 5 ps",\
    '../3_analyze/output/layers_c2_at_ref.dat' i 2 u 1:2:3 w yerrorbars ls 10023 notitle

#====
#plot
#====
set origin 1.05, 0 
set size 1.1, 1.4
set xrange [xmin : xmax]
#set yrange [1.3 : 2.3]
#set yrange [2.0 : 4.5]
unset yrange
set ytics 0.2 
set xlabel "d (\305)" 
set border 1+2+4+8
set label 1 "(b)" left at graph 0.01, graph 0.95
set ylabel "{/Symbol t}_2 (ps)"
plot \
      '../3_analyze/output/tau2.dat' u 1:2 with line ls 1000 notitle,\
      '../3_analyze/output/tau2.dat' u 1:2:3 with yerrorbars ls 1000 notitle

unset multiplot
set term wxt
