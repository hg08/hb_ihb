load 'colorstyle.gnuplot'
load '../2_orientation/output/SYSTEM_orientation_info.txt
set term postscript eps color solid linewidth 2 "Arial" 36 enhanced
set output output_plot_fit_para 
set encoding iso_8859_1

set size square 3.10,1.4

set multiplot
# Plot c(t*) of d
xmin=0.95
xmax=6.05
shift=1.0


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
yc1=0.4
yc2=0.27
yc5=0.10
set arrow 1 from xmin, yc1 to xmax, yc1 nohead ls 20021
set arrow 2 from xmin, yc2 to xmax, yc2 nohead ls 20022
set arrow 3 from xmin, yc5 to xmax, yc5 nohead ls 20023
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
set origin shift, 0 
set size 1.1, 1.4
set xrange [xmin : xmax]
set yrange [0 : 3.0]
set ytics 0.5
set xlabel "d (\305)" 
set border 1+2+4+8
set label 1 "(b)" left at graph 0.01, graph 0.95
set ylabel "{/Symbol t}_i (ps)"
set key at graph 0.85, graph 0.55
tc1=0.32
tc2=2.31
unset arrow 3
set arrow 1 from xmin, tc1 to xmax, tc1 nohead ls 20021
set arrow 2 from xmin, tc2 to xmax, tc2 nohead ls 20022

plot \
      input_plot_fit_para every 4::1 u ($0+1):2 with line ls 10021 title "{/Symbol t}_1",\
      input_plot_fit_para every 4::1 u ($0+1):2:4 with yerrorbars ls 10021 notitle,\
      input_plot_fit_para every 4::3 u ($0+1):2 with line ls 10022 title "{/Symbol t}_{2}",\
      input_plot_fit_para every 4::3 u ($0+1):2:4 with yerrorbars ls 10022 notitle
#====
#plot
#====

set origin 2*shift, 0 
set size 1.1, 1.4
set xrange [xmin : xmax]
set yrange [0 : 0.8]
set ytics 0.2
set xlabel "d (\305)" 
set border 1+2+4+8
set label 1 "(c)" left at graph 0.01, graph 0.95
set ylabel "A_i"
set key at graph 0.85, graph 0.55
ac1 = 0.12
ac2 = 0.61
set arrow 1 from xmin, ac1 to xmax, ac1 nohead ls 20021
set arrow 2 from xmin, ac2 to xmax, ac2 nohead ls 20022
plot \
      input_plot_fit_para every 4::0 u ($0+1):2 with line ls 10021 title "A_1",\
      input_plot_fit_para every 4::0 u ($0+1):2:4 with yerrorbars ls 10021 notitle,\
      input_plot_fit_para every 4::2 u ($0+1):2 with line ls 10022 title "A_2",\
      input_plot_fit_para every 4::2 u ($0+1):2:4 with yerrorbars ls 10022 notitle

unset multiplot
set term wxt