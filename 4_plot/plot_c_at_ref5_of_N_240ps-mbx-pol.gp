load 'colorstyle.gnuplot'

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
set output "c_at_ref5_of_N.eps"
set encoding iso_8859_1

set size square 2.10,1.4

set multiplot
# Plot c(t*) of d
xmin=122
xmax=515

set xrange [xmin : xmax]
set yrange [0 : 0.6]
set xtics 100 
set mxtics 2
set ytics 0.20 
set mytics 2
set format y "%2.1f"
#====
#plot A part 1
#====
set origin 0, 0 
set size 1.1, 1.4
set xlabel 'Number of water molecules' 
set border 1+2+4+8
#set ylabel "{/Symbol t}@^{(s)}_a (ps)"
set xtics (125,216,343,512)
set xlabel "Number of water molecules"
set ylabel "c(t^*=5)"
#set label 2 "ADH criterion" right at graph 0.40, graph 0.95
set label 1 "(c)" left at graph 0.01, graph 0.95
set label 2 "scenario 1 (LC)" right at graph 0.95, graph 0.95
set key left at graph 0.07, graph 0.9
#set arrow 1 from 2.8, graph 0.95 to 3.6, graph 0.95 nohead ls 20000 
#set arrow 2 from 2.8, graph 0.85 to 3.6, graph 0.85 nohead ls 10000
plot [xmin:xmax] \
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_1A.dat" u 1:3 with line ls 10000 title "1 \305",\
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_1A.dat" u 1:3:4 with yerrorbars ls 10000 notitle,\
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_2A.dat" u 1:3 with line ls 10021 title "2 \305",\
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_2A.dat" u 1:3:4 with yerrorbars ls 10021 notitle,\
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_3A.dat" u 1:3 with line ls 10023 title "3 \305",\
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_3A.dat" u 1:3:4 with yerrorbars ls 10023 notitle
#====
#plot A part 2
#====
set origin 0, 0 
set size 1.1, 1.4
set xlabel 'Number of water molecules' 
set border 1+2+4+8
#set ylabel "{/Symbol t}@^{(s)}_a (ps)"
set xtics (125,216,343,512)
set xlabel "Number of water molecules"
#set ylabel "c(t^*=5)"
#set label 2 "ADH criterion" right at graph 0.40, graph 0.95
#set label 1 "(a)" left at graph 0.01, graph 0.95
set label 2 "scenario 1 (LC)" right at graph 0.95, graph 0.95
set key left at graph 0.40, graph 0.9
#set arrow 1 from 2.8, graph 0.95 to 3.6, graph 0.95 nohead ls 20000 
#set arrow 2 from 2.8, graph 0.85 to 3.6, graph 0.85 nohead ls 10000
plot [xmin:xmax] \
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_4A.dat" u 1:3 with line ls 10027 title "4 \305",\
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_4A.dat" u 1:3:4 with yerrorbars ls 10027 notitle,\
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_5A.dat" u 1:3 with line ls 10028 title "5 \305",\
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_5A.dat" u 1:3:4 with yerrorbars ls 10028 notitle,\
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_6A.dat" u 1:3 with line ls 10029 title "6 \305",\
"../3_analyze/output/240ps-mbx-pol_1_case1_c_at_ref5_6A.dat" u 1:3:4 with yerrorbars ls 10029 notitle

#======
#plot B
#======
unset key
shift = 1.0
set origin shift, 0 
set size 1.1, 1.4
set xlabel 'Number of water molecules' 
set border 1+2+4+8
set ylabel "c^{(s)}(t^*=5)"
set label 1 "(d)" left at graph 0.01, graph 0.95
set label 2 "scenario 2 (IHB)" right at graph 0.95, graph 0.95
#set arrow 1 from 2.8, graph 0.95 to 3.6, graph 0.95 nohead ls 20000 
#set arrow 2 from 2.8, graph 0.85 to 3.6, graph 0.85 nohead ls 10000
plot [xmin:xmax] \
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_1A.dat" u 1:3 with line ls 10000 title "1 \305",\
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_1A.dat" u 1:3:4 with yerrorbars ls 10000 notitle,\
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_2A.dat" u 1:3 with line ls 10021 title "2 \305",\
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_2A.dat" u 1:3:4 with yerrorbars ls 10021 notitle,\
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_3A.dat" u 1:3 with line ls 10023 title "3 \305",\
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_3A.dat" u 1:3:4 with yerrorbars ls 10023 notitle,\
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_4A.dat" u 1:3 with line ls 10027 title "4 \305",\
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_4A.dat" u 1:3:4 with yerrorbars ls 10027 notitle,\
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_5A.dat" u 1:3 with line ls 10028 title "5 \305",\
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_5A.dat" u 1:3:4 with yerrorbars ls 10028 notitle,\
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_6A.dat" u 1:3 with line ls 10029 title "6 \305",\
"../3_analyze/output/240ps-mbx-pol_2_case2_c_at_ref5_6A.dat" u 1:3:4 with yerrorbars ls 10029 notitle

unset multiplot
set term wxt
