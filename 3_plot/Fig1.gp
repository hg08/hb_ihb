#load '/home/dylan/bin/vdos_template.gnuplot'
load 'colorstyle.gnuplot'
set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output "Fig1.eps"
set size square 2.2,1.4

set multiplot layout 2,2

xmax=5

set xrange [0 : xmax]
set yrange [0.0: 1.0] 
set xtics 1 
set ytics 0.2 
set mytics 2
set format y "%3.1f"

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
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_1.0.dat' u 1:2:3 w yerrorbars  ls 10021 notitle,\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_1.0.dat' u 1:2 w l ls 10021 title "1 \305",\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_2.0.dat' u 1:2:3 w yerrorbars  ls 10022 notitle,\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_2.0.dat' u 1:2 w l ls 10022 title "2 \305",\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_3.0.dat' u 1:2:3 w yerrorbars  ls 10023 notitle,\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_3.0.dat' u 1:2 w l ls 10023 title "3 \305",\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_4.0.dat' u 1:2:3 w yerrorbars  ls 10024 notitle,\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_4.0.dat' u 1:2 w l ls 10024 title "4 \305",\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_5.0.dat' u 1:2:3 w yerrorbars  ls 10028 notitle,\
     '../1_case1/output/128w_itp_wat_pair_hbacf_h_ihb_5.0.dat' u 1:2 w l ls 10028 title "5 \305"

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
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_1.dat' u 1:2:3 w yerrorbars  ls 10021 notitle,\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_1.dat' u 1:2 w l ls 10021 title "1 \305",\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_2.dat' u 1:2:3 w yerrorbars  ls 10022 notitle,\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_2.dat' u 1:2 w l ls 10022 title "2 \305",\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_3.dat' u 1:2:3 w yerrorbars  ls 10023 notitle,\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_3.dat' u 1:2 w l ls 10023 title "3 \305",\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_4.dat' u 1:2:3 w yerrorbars  ls 10024 notitle,\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_4.dat' u 1:2 w l ls 10024 title "4 \305",\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_5.dat' u 1:2:3 w yerrorbars  ls 10028 notitle,\
     '../2_case2/output/128w_itp_wat_pair_hbacf_h_ihb_5.dat' u 1:2 w l ls 10028 title "5 \305"


unset multiplot
set term wxt
