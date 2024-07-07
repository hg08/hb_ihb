load 'colorstyle.gnuplot'
set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output "Fig2.eps"
set encoding iso_8859_1
set size square 2.2,1.4

set multiplot layout 2,2

xmax=5

set xrange [0 : xmax]
set yrange [0.0: 1.0] 
set xtics 1 
set ytics 0.2 
set mytics 2
set format y "%3.1f"

# Add arrows and labels for both plots
dash_length = 4  # Length of the dash
dash_space = 1    # Space between dashes, adjust this to change the spacing
set arrow 1 from 2, graph 0.2 to 2, graph 0.8 nohead dashtype (dash_length, dash_space) linecolor rgb "#878787" linewidth 5 front

# Set up the triangle as an arrowhead
triangle_width = 0.08  # Adjust this value to change the width of the triangle
set object 1 polygon from \
     2, graph 0.82 to \
     2 + triangle_width, graph 0.78 to \
     2 - triangle_width, graph 0.78 to \
     2, graph 0.82
set object 1 fc rgb "#878787" fillstyle solid noborder front
set label 3 at 2, graph 0.17 "increase d" center font "Arial,32" textcolor rgb "black"

plot_interval = 3  # Change this value to adjust the plotting frequency

#====
#plot A
# ADH scheme 1
#====
set origin 0, 0.0
set size 1.1, 1.4
set border 1+2+4+8
set ylabel "c(t)"
set xlabel "t (ps)"
set label 1 "(a)" right at graph 0.1, graph 0.95
set label 2 "scenario 1 (LC)" left at graph 0.13, graph 0.95
#set arrow 1 ls 4 from graph 0.4, graph 0.25 to graph 0.4, graph 0.75 
plot [0.0:xmax] \
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_1.0.dat' every plot_interval::0 u 1:2 w l ls 10021 title "1 \305",\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_2.0.dat' every plot_interval::0 u 1:2 w l ls 10022 title "2 \305",\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_3.0.dat' every plot_interval::0 u 1:2 w l ls 10023 title "3 \305",\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_4.0.dat' every plot_interval::0 u 1:2 w l ls 10024 title "4 \305",\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_5.0.dat' every plot_interval::0 u 1:2 w l ls 10028 title "5 \305"

#====
#plot B
# ADH scheme2
#====
set origin 1.05, 0.0
set size 1.1, 1.4
set border 1+2+4+8
set xlabel "t (ps)"
set ylabel "c^{(s)}(t)"
set label 1 "(b)" right at graph 0.1, graph 0.95
set label 2 "scenario 2 (IHB)" left at graph 0.13, graph 0.95
plot [0.0:xmax] \
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_1.dat' every plot_interval::0 u 1:2 w l ls 10021 title "1 \305",\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_2.dat' every plot_interval::0 u 1:2 w l ls 10022 title "2 \305",\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_3.dat' every plot_interval::0 u 1:2 w l ls 10023 title "3 \305",\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_4.dat' every plot_interval::0 u 1:2 w l ls 10024 title "4 \305",\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_5.dat' every plot_interval::0 u 1:2 w l ls 10028 title "5 \305"


unset arrow 1
unset object 1
unset label 3

unset multiplot
set term wxt
