load 'colorstyle.gnuplot'
set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output "FigS2.eps"
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
xmin=0.00
xmax=31.2808

set origin 0, 0 
set xlabel "z (\305)" 
set xrange [xmin : xmax]
set xtics 5 
set ytics 0.5
set mytics 2
set format y "%3.1f"

#====
#plot
#====
set size 1.1, 1.4
set border 1+2+4+8
set ylabel "{/Symbol r}(z)  (g/cm^{3})"
#set label 1 "b" right at graph 0.05, graph 0.95
set yrange [-0.01: 1.53] 
#set arrow 3 from 1.0,0.142 to 6,0.142 nohead ls 4 
#set label 7 "0.10" at 0.75,0.25 textcolor rgb '#000000'
plot [xmin:xmax] \
     '../2_density/output/sorted_density_OH.dat' u 1:2 with p ls 1 notitle,\
     '../2_density/output/sorted_density_OH.dat' u 1:2 with l ls 1 notitle

unset multiplot
