load 'colorstyle.gnuplot'
load '../2_density/output/216h2o-240ps-mbx-pol_density_info.txt'

set term postscript eps color solid linewidth 2 "Arial" 34 enhanced
set output "density_125-512.eps"
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


set size square 2.90,1.3

set multiplot
xmin=20.00
xmax=60
shift=0.9

set origin 0, 0 
set xlabel "z (\305)" 
set xrange [xmin : xmax]
set xtics 10 
set ytics 0.5
set mytics 2
set format y "%3.1f"

#====
#plot
#====
set size 1.0, 1.3
set border 1+2+4+8
set ylabel "{/Symbol r}(z)  (g/cm^{3})"
set label 1 "(a)" right at graph 0.1, graph 0.95
set yrange [-0.01: 1.5] 
set label 2 "125 water molecules" at graph 0.15, graph 0.95 textcolor rgb '#000000'

plot [xmin:xmax] \
     '../2_density/output/density_OH_125h2o-240ps-mbx-pol.dat' u 1:2 with p ls 1 notitle,\
     '../2_density/output/density_OH_125h2o-240ps-mbx-pol.dat' u 1:2 with l ls 1 notitle

#====
#plot
#====
xmin=30.00
xmax=70
set origin shift, 0 
set size 1.0, 1.3
set border 1+2+4+8
unset ylabel 
set label 1 "(b)" right at graph 0.1, graph 0.95
set xrange [xmin : xmax]
set yrange [-0.01: 1.5] 
#set arrow 1 from xmin,1.00 to xmax,1.00 nohead ls 4 
set label 2 "216 water molecules" at graph 0.15, graph 0.95 textcolor rgb '#000000'

plot [xmin:xmax] \
     input_plot u 1:2 with p ls 1 notitle,\
     input_plot u 1:2 with l ls 1 notitle

#====
#plot
#====
xmin=30
xmax=90
set origin 2*shift, 0 
set size 1.0, 1.3
set border 1+2+4+8
unset ylabel 
set label 1 "(c)" right at graph 0.1, graph 0.95
set xrange [xmin : xmax]
set yrange [-0.01: 1.5] 
#set arrow 1 from xmin,1.00 to xmax,1.00 nohead ls 4 
set label 2 "512 water molecules" at graph 0.15, graph 0.95 textcolor rgb '#000000'

plot [xmin:xmax] \
     '../2_density/output/density_OH_512h2o-240ps-mbx-pol.dat' u 1:2 with p ls 1 notitl,\
     '../2_density/output/density_OH_512h2o-240ps-mbx-pol.dat' u 1:2 with l ls 1 notitle

unset multiplot
