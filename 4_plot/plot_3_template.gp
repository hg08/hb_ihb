load 'colorstyle.gnuplot'
load '../2_orientation/output/SYSTEM_orientation_info.txt

set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output output_cHB_ref
set encoding iso_8859_1

set size square 1.10,1.4

set multiplot
# Plot c(t*) of d
xmin=0.95
xmax=6.05

set xrange [xmin : xmax]
set yrange [0 : 1.0]
set xtics 1 
#set mxtics 2
set ytics 0.20 
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
set ylabel "{c}(t*), {c}^{(s)}(t*)"
#set label 2 "ADH criterion" right at graph 0.40, graph 0.95
plot [xmin:xmax] \
             input_case1_cHB_ref i 0 u 1:2 with line ls 20021 title "scenario 1: t* = 1 ps",\
             input_case1_cHB_ref i 0 u 1:2:3 with yerrorbars ls 10021 notitle,\
             input_case1_cHB_ref i 1 u 1:2 with line ls 20022 title "                 2 ps",\
             input_case1_cHB_ref i 1 u 1:2:3 with yerrorbars ls 10022 notitle,\
             input_case1_cHB_ref i 2 u 1:2 with line ls 20023 title "                 5 ps",\
             input_case1_cHB_ref i 2 u 1:2:3 with yerrorbars ls 10023 notitle,\
             input_case2_cHB_ref i 0 u 1:2 with line ls 10021 title "scenario 2: t*= 1 ps",\
             input_case2_cHB_ref i 0 u 1:2:3 with yerrorbars ls 10021 notitle,\
             input_case2_cHB_ref i 1 u 1:2 with line ls 10022 title "                2 ps",\
             input_case2_cHB_ref i 1 u 1:2:3 with yerrorbars ls 10022 notitle,\
             input_case2_cHB_ref i 2 u 1:2 with line ls 10023 title "                5 ps",\
             input_case2_cHB_ref i 2 u 1:2:3 with yerrorbars ls 10023 notitle
unset multiplot
set term wxt
