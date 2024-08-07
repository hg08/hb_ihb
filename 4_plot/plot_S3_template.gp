load 'colorstyle.gnuplot'
load '../2_orientation/output/SYSTEM_orientation_info.txt

set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output output_plot  # c2 without fitting
set encoding iso_8859_1

set size square 1.10,1.4

set multiplot
xmin=0.99
xmax=6.01

set xrange [xmin : xmax]
set mytics 2
set format y "%3.1f"

set xrange [0 : 5]
set xtics 1
set ytics 0.2
#====
#plot
#====
set origin 0, 0
set size 1.1, 1.4
set xlabel "t (ps)"
set border 1+2+4+8
set ylabel "C_2(t)"
#set label 1 "a" right at graph 0.05, graph 0.95
set yrange [-0.: 1.0]
plot input_plot1 u 1:2 w l ls 10021 title "1 \305",\
     input_plot2 u 1:2 w l ls 10022 title "2 \305",\
     input_plot3 u 1:2 w l ls 10023 title "3 \305",\
     input_plot4 u 1:2 w l ls 10024 title "4 \305",\
     input_plot5 u 1:2 w l ls 10025 title "5 \305",\
     input_plot6 u 1:2 w l ls 10026 title "6 \305"

unset multiplot
