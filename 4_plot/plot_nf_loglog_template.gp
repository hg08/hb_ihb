load 'colorstyle.gnuplot'
load '../2_orientation/output/SYSTEM_orientation_info.txt'

set style line 30021 lt 1 lc rgb '#5e4fa2' lw 2 pt 1 ps 1
set style line 30022 lt 1 lc rgb '#3288bd' lw 2 pt 1 ps 1
set style line 30023 lt 1 lc rgb '#66c2a5' lw 2 pt 1 ps 1
set style line 30024 lt 1 lc rgb '#abdaa4' lw 2 pt 8 ps 1
set style line 30025 lt 1 lc rgb '#e6f598' lw 2 pt 8 ps 1
set style line 30026 lt 1 lc rgb '#fee08b' lw 2 pt 7 ps 1
set style line 30027 lt 1 lc rgb '#fdae61' lw 2 pt 12 ps 1
set style line 30028 lt 1 lc rgb '#f46d43' lw 2 pt 13 ps 1
set style line 30029 lt 1 lc rgb '#d53e4f' lw 2 pt 4 ps 1
set style line 30030 lt 1 lc rgb '#9e0142' lw 2 pt 8 ps 1

set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output output_nf_loglog 
set encoding iso_8859_1
set size square 2.2,1.4

set multiplot layout 2,2

xmin=0.03
xmax=10
ymin=0.02
ymax=1.1

set xrange [xmin : xmax]
set yrange [ymin : ymax] 
#set xtics 1 
#set ytics 0.2 
#set mytics 2
set logscale xy
set format y "%3.1f"

## Add arrows and labels for both plots
#dash_length = 4  # Length of the dash
#dash_space = 1    # Space between dashes, adjust this to change the spacing
#set arrow 1 from 2, graph 0.1 to 2, graph 0.65 nohead dashtype (dash_length, dash_space) linecolor rgb "#878787" linewidth 5 front

## Set up the triangle as an arrowhead
#triangle_width = 0.08  # Adjust this value to change the width of the triangle
#y_coord1=0.67
#y_coord2=0.63
#set object 1 polygon from \
#     2, graph y_coord1 to \
#     2 + triangle_width, graph y_coord2 to \
#     2 - triangle_width, graph y_coord2 to \
#     2, graph y_coord1
#set object 1 fc rgb "#878787" fillstyle solid noborder front
#set label 3 at 2, graph 0.06 "increase d" center font "Arial,32" textcolor rgb "black"

plot_interval = 3  # Change this value to adjust the plotting frequency


#====
#plot A
# ADH scheme 1
#====
set origin 0, 0.0
set size 1.1, 1.4
set border 1+2+4+8
set ylabel "n_@f^{(s)}(t)"
set xlabel "t (ps)"
set label 1 "(a)" right at graph 0.99, graph 0.95
set label 2 "R-{/Symbol b} def." left at graph 0.03, graph 0.65
set key left at graph 0.02, graph 0.55
#set arrow 1 ls 4 from graph 0.4, graph 0.25 to graph 0.4, graph 0.75 
# The variable system is passed from the command line
exp = 2.7182818
k1 = -0.27
c1 = 0.368
f1(x) = c1 * (x**k1)
k2 = -1.19
c2 = 0.7
f2(x) = c2 * (x**k2)
k3 = -0.13
c3 = 0.60
f3(x) = c3 * (x**k3)
k4 = -1.19
c4 = 2.3
f4(x) = c4 * (x**k4)
k5 = -0.10
c5 = 0.72
f5(x) = c5 * (x**k5)
k6 = -1.19
c6 = 3.3
f6(x) = c6 * (x**k6)
plot [xmin:xmax] \
     input_case2_nf_1_1 every plot_interval::1 u 1:2 w lp ls 10021 notitle "1 \305",\
     input_case2_nf_1_2 every plot_interval::1 u 1:2 w lp ls 10022 notitle "2 \305",\
     input_case2_nf_1_3 every plot_interval::1 u 1:2 w lp ls 10023 notitle "3 \305",\
     input_case2_nf_1_4 every plot_interval::1 u 1:2 w lp ls 10024 notitle "4 \305",\
     input_case2_nf_1_5 every plot_interval::1 u 1:2 w lp ls 10028 notitle "5 \305",\
     input_case2_nf_1_6 every plot_interval::1 u 1:2 w lp ls 10029 notitle "6 \305",\
     input_case2_nf_1_1 every plot_interval::1 u 1:2:3 w errorbars ls 10021 title "1 \305",\
     input_case2_nf_1_2 every plot_interval::1 u 1:2:3 w errorbars ls 10022 title "2 \305",\
     input_case2_nf_1_3 every plot_interval::1 u 1:2:3 w errorbars ls 10023 title "3 \305",\
     input_case2_nf_1_4 every plot_interval::1 u 1:2:3 w errorbars ls 10024 title "4 \305",\
     input_case2_nf_1_5 every plot_interval::1 u 1:2:3 w errorbars ls 10028 title "5 \305",\
     input_case2_nf_1_6 every plot_interval::1 u 1:2:3 w errorbars ls 10029 title "6 \305"

#====
#plot B
# ADH scheme2
#====
set origin 1.05, 0.0
set size 1.1, 1.4
set border 1+2+4+8
set xlabel "t (ps)"
set ylabel "n_@f^{(s)}(t)"
set label 1 "(b)" right at graph 0.99, graph 0.95
set label 2 "R-{/Symbol q} def." left at graph 0.03, graph 0.65
k1 = -0.333
c1 = 0.31
f1(x) = c1 * (x**k1)
k2 = -1.19
c2 = 0.45
f2(x) = c2 * (x**k2)
k3 = -0.18
c3 = 0.55
f3(x) = c3 * (x**k3)
k4 = -1.19
c4 = 1.5
f4(x) = c4 * (x**k4)
k5 = -0.12
c5 = 0.70
f5(x) = c5 * (x**k5)
k6 = -1.19
c6 = 3.0
f6(x) = c6 * (x**k6)
plot [xmin:xmax] \
     input_case2_nf_2_1 every plot_interval::1 u 1:2 w lp ls 10021 notitle "1 \305",\
     input_case2_nf_2_2 every plot_interval::1 u 1:2 w lp ls 10022 notitle "2 \305",\
     input_case2_nf_2_3 every plot_interval::1 u 1:2 w lp ls 10023 notitle "3 \305",\
     input_case2_nf_2_4 every plot_interval::1 u 1:2 w lp ls 10024 notitle "4 \305",\
     input_case2_nf_2_5 every plot_interval::1 u 1:2 w lp ls 10028 notitle "5 \305",\
     input_case2_nf_2_6 every plot_interval::1 u 1:2 w lp ls 10029 notitle,\
     input_case2_nf_2_1 every plot_interval::1 u 1:2:3 w errorbars ls 10021 title "1 \305",\
     input_case2_nf_2_2 every plot_interval::1 u 1:2:3 w errorbars ls 10022 title "2 \305",\
     input_case2_nf_2_3 every plot_interval::1 u 1:2:3 w errorbars ls 10023 title "3 \305",\
     input_case2_nf_2_4 every plot_interval::1 u 1:2:3 w errorbars ls 10024 title "4 \305",\
     input_case2_nf_2_5 every plot_interval::1 u 1:2:3 w errorbars ls 10028 title "5 \305",\
     input_case2_nf_2_6 every plot_interval::1 u 1:2:3 w errorbars ls 10029 title "6 \305"
     #[0.04: 2] f1(x) w l ls 4 notitle,\
     #[1.3: 10] f2(x) w l ls 20029 notitle,\
     #[0.04: 3.4] f3(x) w l ls 4 notitle,\
     #[2.2: 10] f4(x) w l ls 20029 notitle,\
     #[0.04: 4.5] f5(x) w l ls 4 notitle,\
     #[3.2: 10] f6(x) w l ls 20029 notitle


unset arrow 1
unset object 1
unset label 3

unset multiplot
set term wxt
