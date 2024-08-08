load 'colorstyle.gnuplot'
load '../2_orientation/output/512h2o-240ps-mbx-pol_orientation_info.txt

#========================
#渐变2A：blue-red-yellow
#========================
set style line 10021 lt 1 lc rgb '#5e4fa2' lw 2
set style line 10022 lt 1 lc rgb '#3288bd' lw 2
set style line 10023 lt 1 lc rgb '#66c2a5' lw 2
set style line 10024 lt 1 lc rgb '#abdaa4' lw 2
set style line 10025 lt 1 lc rgb '#e6f598' lw 2
set style line 10026 lt 1 lc rgb '#fee08b' lw 2
set style line 10027 lt 1 lc rgb '#fdae61' lw 2
set style line 10028 lt 1 lc rgb '#f46d43' lw 2
set style line 10029 lt 1 lc rgb '#d53e4f' lw 2
set style line 10030 lt 1 lc rgb '#9e0142' lw 2
#========================
#渐变2B：blue-red-yellow
#========================
set style line 20021 linecolor rgb '#5e4fa2' lt 1 linewidth 2 dashtype '.-'
set style line 20022 linecolor rgb '#3288bd' lt 1 linewidth 2 dashtype '.-'
set style line 20023 linecolor rgb '#66c2a5' lt 1 linewidth 2 dashtype '.-'
set style line 20024 linecolor rgb '#abdaa4' lt 1 linewidth 2 dashtype '.-'
set style line 20025 linecolor rgb '#e6f598' lt 1 linewidth 2 dashtype '.-'
set style line 20026 linecolor rgb '#fee08b' lt 1 linewidth 2 dashtype '.-'
set style line 20027 linecolor rgb '#fdae61' lt 1 linewidth 2 dashtype '.-'
set style line 20028 linecolor rgb '#f46d43' lt 1 linewidth 2 dashtype '.-'
set style line 20029 linecolor rgb '#d53e4f' lt 1 linewidth 2 dashtype '.-'
set style line 20030 linecolor rgb '#9e0142' lt 1 linewidth 2 dashtype '.-'
#========================
#渐变2A：blue-red-yellow
#========================
set style line 30021 lt 1 lc rgb '#5e4fa2' lw 2 pt 1 ps 1
set style line 30022 lt 1 lc rgb '#3288bd' lw 2 pt 1 ps 1
set style line 30023 lt 1 lc rgb '#66c2a5' lw 2 pt 1 ps 1
set style line 30024 lt 1 lc rgb '#abdaa4' lw 1 pt 8 ps 1.5
set style line 30025 lt 1 lc rgb '#e6f598' lw 1 pt 8 ps 1.5
set style line 30026 lt 1 lc rgb '#fee08b' lw 1 pt 7 ps 1.5
set style line 30027 lt 1 lc rgb '#fdae61' lw 1 pt 12 ps 1.5
set style line 30028 lt 1 lc rgb '#f46d43' lw 1 pt 13 ps 1.5
set style line 30029 lt 1 lc rgb '#d53e4f' lw 1 pt 4 ps 1.5
set style line 30030 lt 1 lc rgb '#9e0142' lw 1 pt 8 ps 1.5

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
             input_case1_cHB_ref i 0 u 1:2 with lp ls 20021 notitle,\
             input_case1_cHB_ref i 0 u 1:2:3 with yerrorbars ls 10021 notitle,\
             input_case1_cHB_ref i 1 u 1:2 with lp ls 20022 notitle,\
             input_case1_cHB_ref i 1 u 1:2:3 with yerrorbars ls 10022 notitle,\
             input_case1_cHB_ref i 2 u 1:2 with lp ls 20023 notitle,\
             input_case1_cHB_ref i 2 u 1:2:3 with yerrorbars ls 10023 notitle,\
             input_case2_cHB_ref i 0 u 1:2 with lp ls 30021 title "t*= 1 ps",\
             input_case2_cHB_ref i 0 u 1:2:3 with yerrorbars ls 10021 notitle,\
             input_case2_cHB_ref i 1 u 1:2 with lp ls 30022 title "    2 ps",\
             input_case2_cHB_ref i 1 u 1:2:3 with yerrorbars ls 10022 notitle,\
             input_case2_cHB_ref i 2 u 1:2 with lp ls 30023 title "    5 ps",\
             input_case2_cHB_ref i 2 u 1:2:3 with yerrorbars ls 10023 notitle
unset multiplot
set term wxt
