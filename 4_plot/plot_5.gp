load 'colorstyle.gnuplot'
load '../2_orientation/output/512h2o-240ps-mbx-pol_orientation_info.txt'

set term postscript eps color solid linewidth 2 "Arial" 36 enhanced
set output 'Fig5.eps'
set encoding iso_8859_1

set size square 2.10,2.8

set multiplot
# Plot c(t*) of d
xmin=0.95
xmax=6.05
shift=1.0
yshift=1.4

#====
#plot
#====
set origin 0, yshift 
set size 1.1, 1.4
set xlabel "t (ps)" 
set border 1+2+4+8
set ylabel "C_{2}(t)"
set format y "%2.1f"
set label 1 "(a)" right at graph 0.1, graph 0.95
set xrange [-0.05: 5.05] 
set yrange [0.0: 1.0] 


f1(x)=a1*exp(-x/b1)+c1*exp(-x/d1)
f2(x)=a2*exp(-x/b2)+c2*exp(-x/d2)
f3(x)=a3*exp(-x/b3)+c3*exp(-x/d3)
f4(x)=a4*exp(-x/b4)+c4*exp(-x/d4)
f5(x)=a5*exp(-x/b5)+c5*exp(-x/d5)
f6(x)=a6*exp(-x/b6)+c6*exp(-x/d6)

fit [0:5] f1(x) input_plot1 u 1:2 via a1,b1,c1,d1 
fit [0:5] f2(x) input_plot2 u 1:2 via a2,b2,c2,d2 
fit [0:5] f3(x) input_plot3 u 1:2 via a3,b3,c3,d3 
fit [0:5] f4(x) input_plot4 u 1:2 via a4,b4,c4,d4 
fit [0:5] f5(x) input_plot5 u 1:2 via a5,b5,c5,d5 
fit [0:5] f6(x) input_plot6 u 1:2 via a6,b6,c6,d6 

plot input_plot1 every 3::1 u 1:2 w p ls 30021 title "1 \305",\
 input_plot4 every 3::1 u 1:2 w p ls 30023 title "4 \305",\
 input_plot6 every 3::1 u 1:2 w p ls 30029 title "6 \305",\
 f1(x) w l ls 10021 title "1 \305 (fitted)",\
 f4(x) w l ls 10023 title "4 \305 (fitted)",\
 f6(x) w l ls 10029 title "6 \305 (fitted)"

#====
#plot
#====
set origin shift, yshift 
set size 1.1, 1.4
set xlabel "d (\305)" 
set border 1+2+4+8
set xrange [xmin : xmax]
set yrange [0 : 0.6]
set xtics 1 
set ytics 0.2 
set ylabel "C_2(t*)"
set key at graph 0.8, graph 0.9
set label 1 "(b)" left at graph 0.01, graph 0.95
yc1=0.385
yc2=0.240
yc5=0.077
set arrow 1 from xmin, yc1 to xmax, yc1 nohead ls 20021
set arrow 2 from xmin, yc2 to xmax, yc2 nohead ls 20022
set arrow 3 from xmin, yc5 to xmax, yc5 nohead ls 20023
plot \
    input_c2_ref i 0 u 1:2 w l ls 10021 title "t* = 1 ps",\
    input_c2_ref i 0 u 1:2:3 w yerrorbars ls 10021 notitle,\
    input_c2_ref i 1 u 1:2 w l ls 10022 title "t* = 2 ps",\
    input_c2_ref i 1 u 1:2:3 w yerrorbars ls 10022 notitle,\
    input_c2_ref i 2 u 1:2 w l ls 10023 title "t* = 5 ps",\
    input_c2_ref i 2 u 1:2:3 w yerrorbars ls 10023 notitle

#====
#plot
#====
set origin 0, 0 
set size 1.1, 1.4
set xrange [xmin : xmax]
set yrange [0 : 3.0]
set ytics 0.5
set xlabel "d (\305)" 
set border 1+2+4+8
set label 1 "(c)" left at graph 0.01, graph 0.95
set ylabel "{/Symbol t}_i (ps)"
set key at graph 0.85, graph 0.55
tc1=0.32
tc2=2.22
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

set origin shift, 0 
set size 1.1, 1.4
set xrange [xmin : xmax]
set yrange [0 : 0.8]
set ytics 0.2
set xlabel "d (\305)" 
set border 1+2+4+8
set label 1 "(d)" left at graph 0.01, graph 0.95
set ylabel "A_i"
set key at graph 0.85, graph 0.55
ac1 = 0.135
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
