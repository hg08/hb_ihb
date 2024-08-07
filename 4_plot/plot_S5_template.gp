load 'colorstyle.gnuplot'
load '../2_orientation/output/SYSTEM_orientation_info.txt
set encoding iso_8859_1 #

set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output output_plot_single_fit
set size square 1.10,1.4

set multiplot
set key spacing 1.2

set xrange [0 : 5]
set xtics 1 
#set mxtics 2
set ytics 0.2 
set mytics 2
#set format y "%.1t"
#set format y "%3.2f
set format y "%2.1f
#====
#plot
#====
set origin 0, 0 
set size 1.1, 1.4
set xlabel "t (ps)" 
set border 1+2+4+8
set ylabel "C_{2}(t)"
#set label 1 "a" right at graph 0.05, graph 0.95
set yrange [-0.: 1.0] 


f(x)=exp(-x/b1)
f2(x)=exp(-x/b2)
f3(x)=exp(-x/b3)
f4(x)=exp(-x/b4)
f5(x)=exp(-x/b5)
f6(x)=exp(-x/b6)

fit [0.2:5] f(x) input_plot1 u 1:2 via b1 
fit [0.2:5] f2(x) input_plot2 u 1:2 via b2 
fit [0.2:5] f3(x) input_plot3 u 1:2 via b3 
fit [0.2:5] f4(x) input_plot4 u 1:2 via b4 
fit [0.2:5] f5(x) input_plot5 u 1:2 via b5 
fit [0.2:5] f6(x) input_plot6 u 1:2 via b6 

plot input_plot1 u 1:2 w l ls 10021 title "1 \305",\
 input_plot2 u 1:2 w l ls 10022 title "2 \305",\
 input_plot3 u 1:2 w l ls 10023 title "3 \305",\
 input_plot4 u 1:2 w l ls 10024 title "4 \305",\
 input_plot5 u 1:2 w l ls 10025 title "5 \305",\
 input_plot6 u 1:2 w l ls 10026 title "6 \305",\
 f(x) w l ls 20021 notitle,\
 f2(x) w l ls 20022 notitle,\
 f3(x) w l ls 20023 notitle,\
 f4(x) w l ls 20024 notitle,\
 f5(x) w l ls 20025 notitle,\
 f6(x) w l ls 20026 notitle

