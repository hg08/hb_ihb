load 'colorstyle.gnuplot'
load '../2_orientation/output/128w-120ps-dt40_orientation_info.txt

set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output output_plot
set encoding iso_8859_1

set style line 1000 lt 1 lw 2 lc rgb '#000000' # Black solid
set style line 1 lt 1 lw 1 lc rgb '#636363' # solid
set style line 2 lt 1 lw 1 lc rgb '#bdbdbd' # solid
set style line 3 lt 1 lw 1 lc rgb '#f0f0f0' # solid

set style line 4 lw 1 lc rgb "#000000" lt 3 dashtype '-'
set style line 5 lw 1 lc rgb "#636363" lt 3 dashtype '-'
set style line 6 lw 1 lc rgb "#bdbdbd" lt 3 dashtype '-'
set style line 7 lw 1 lc rgb "#f0f0f0" lt 3 dashtype '-'

set style line 8 lw 2 lc rgb "#000000" lt 3 dashtype '.-'
set style line 9 lw 2 lc rgb "#636363" lt 3 dashtype '.-'
set style line 10 lw 2 lc rgb "#bdbdbd" lt 3 dashtype '.-'
set style line 11 lw 2 lc rgb "#f0f0f0" lt 3 dashtype '.-'

set style line 12 lw 2 lc rgb "#000000" lt 3 dashtype '..-'
set style line 13 lw 2 lc rgb "#636363" lt 3 dashtype '..-'
set style line 14 lw 2 lc rgb "#bdbdbd" lt 3 dashtype '..-'
set style line 15 lw 2 lc rgb "#f0f0f0" lt 3 dashtype '..-'


set size square 1.10,1.4

set multiplot
xmin=0.99
xmax=6.01

set xrange [xmin : xmax]
set mytics 2
set format y "%3.1f"


set xrange [0 : 5]
set xtics 1
set ytics 0.2
#====
#plot
#====
set origin 0, 0
set size 1.1, 1.4
set xlabel "t (ps)"
set border 1+2+4+8
set ylabel "C_2(t)"
#set label 1 "a" right at graph 0.05, graph 0.95
set yrange [-0.: 1.0]
plot '../results/128w_itp_c2_ihb_ave_1.dat' u 1:2 w l ls 1 title "1 \305",\
 '../results/128w_itp_c2_ihb_ave_2.dat' u 1:2 w l ls 2 title "2 \305",\
 '../results/128w_itp_c2_ihb_ave_3.dat' u 1:2 w l ls 3 title "3 \305",\
 '../results/128w_itp_c2_ihb_ave_4.dat' u 1:2 w l ls 4 title "4 \305",\
 '../results/128w_itp_c2_ihb_ave_5.dat' u 1:2 w l ls 5 title "5 \305",\
 '../results/128w_itp_c2_ihb_ave_6.dat' u 1:2 w l ls 6 title "6 \305"

unset multiplot
