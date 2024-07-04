load 'colorstyle.gnuplot'
set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output "Fig3.eps"
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
             '../3_analyze/output/1_case1_layers_c_at_ref.dat' i 0 u 1:2:3 w yerrorbars ls 10021 notitle,\
             '../3_analyze/output/1_case1_layers_c_at_ref.dat' i 0 u 1:2 with linespoints ls 20021 notitle "t*= 1 ps",\
             '../3_analyze/output/1_case1_layers_c_at_ref.dat' i 1 u 1:2:3 w yerrorbars ls 10022 notitle,\
             '../3_analyze/output/1_case1_layers_c_at_ref.dat' i 1 u 1:2 with linespoints ls 20022 notitle "    2 ps",\
             '../3_analyze/output/1_case1_layers_c_at_ref.dat' i 2 u 1:2:3 w yerrorbars ls 10023 notitle,\
             '../3_analyze/output/1_case1_layers_c_at_ref.dat' i 2 u 1:2 with linespoints ls 20023 notitle "    5 ps",\
             '../3_analyze/output/2_case2_layers_c_at_ref.dat' i 0 u 1:2:3 w yerrorbars ls 10021 notitle,\
             '../3_analyze/output/2_case2_layers_c_at_ref.dat' i 0 u 1:2 with linespoints ls 10021 title "t*= 1 ps",\
             '../3_analyze/output/2_case2_layers_c_at_ref.dat' i 1 u 1:2:3 w yerrorbars ls 10022 notitle,\
             '../3_analyze/output/2_case2_layers_c_at_ref.dat' i 1 u 1:2 with linespoints ls 10022 title "t*= 2 ps",\
             '../3_analyze/output/2_case2_layers_c_at_ref.dat' i 2 u 1:2:3 w yerrorbars ls 10023 notitle,\
             '../3_analyze/output/2_case2_layers_c_at_ref.dat' i 2 u 1:2 with linespoints ls 10023 title "t*= 5 ps"
unset multiplot
set term wxt
