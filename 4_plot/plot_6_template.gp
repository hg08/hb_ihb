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
set output output_nf_ref
set encoding iso_8859_1

set size square 1.10,1.4

set multiplot
# Plot c(t*) of d
xmin=0.95
xmax=6.05

set xrange [xmin : xmax]
set yrange [0 : 0.4]
set xtics 1 
#set mxtics 2
set ytics 0.10 
#set mytics 2
set format y "%2.1f"
#====
#plot
#====
set origin 0, 0 
set size 1.1, 1.4
set xlabel "d (\305)" 
set border 1+2+4+8
#set ylabel "{/Symbol t}@^{(s)}_a (ps)"
set ylabel "n_@{f}^{(s)}(t*)"
#set label 2 "ADH criterion" right at graph 0.40, graph 0.95
set label 1 " R-{/Symbol b} def." left at graph 0.05, graph 0.95
set label 2 " R-{/Symbol q} def." left at graph 0.05, graph 0.85
set arrow 1 from 2.6, graph 0.95 to 3.4, graph 0.95 nohead ls 20000 
set arrow 2 from 2.6, graph 0.85 to 3.4, graph 0.85 nohead ls 10000
plot [xmin:xmax] \
             input_case2_nf_criterion1_ref i 0 u 1:2 with line ls 20021 notitle "R-{/Symbol b}: t* = 1 ps",\
             input_case2_nf_criterion1_ref i 0 u 1:2:3 with yerrorbars ls 10021 notitle,\
             input_case2_nf_criterion1_ref i 1 u 1:2 with line ls 20022 notitle,\
             input_case2_nf_criterion1_ref i 1 u 1:2:3 with yerrorbars ls 10022 notitle,\
             input_case2_nf_criterion1_ref i 2 u 1:2 with line ls 20023 notitle,\
             input_case2_nf_criterion1_ref i 2 u 1:2:3 with yerrorbars ls 10023 notitle,\
             input_case2_nf_criterion2_ref i 0 u 1:2 with lp ls 30021 title "t*= 1 ps",\
             input_case2_nf_criterion2_ref i 0 u 1:2:3 with yerrorbars ls 10021 notitle,\
             input_case2_nf_criterion2_ref i 1 u 1:2 with lp ls 30022 title "    2 ps",\
             input_case2_nf_criterion2_ref i 1 u 1:2:3 with yerrorbars ls 10022 notitle,\
             input_case2_nf_criterion2_ref i 2 u 1:2 with lp ls 30023 title "    5 ps",\
             input_case2_nf_criterion2_ref i 2 u 1:2:3 with yerrorbars ls 10023 notitle
unset multiplot
set term wxt
