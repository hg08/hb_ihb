load 'colorstyle.gnuplot'
load '../2_orientation/output/SYSTEM_orientation_info.txt

set term postscript eps color solid linewidth 2 "Arial" 32 enhanced
set output output_kprime
set encoding iso_8859_1

set size square 1.10,1.4

set multiplot
xmin=0.95
xmax=6.05

set origin 0, 0 
set xlabel "d (\305)" 
set xrange [xmin : xmax]
set xtics 1 
set mxtics 2
set ytics 0.1 
set mytics 2
set format y "%3.1f"

#====
#plot
#====
set size 1.1, 1.4
set border 1+2+4+8
set ylabel "k@{/Symbol \242}_{LC}, k@{/Symbol \242}_{IHB}  (ps^{-1})"
set yrange [0.0: 0.3] 
#set arrow 3 from 1.0,0.10 to 6,0.10 nohead ls 4 
#set label 6 "k@{/Symbol \242}_{bulk}=0.10 ps^{-1}" right at graph 0.95, graph 0.40
#set label 7 "0.10" at 0.25,0.10 textcolor rgb '#000000'
plot [xmin:xmax] \
     input_case1_kkprime u 1:4 with line ls 20021 title "scenario 1  (LC)",\
     input_case2_kkprime u 1:4 with line ls 10022 title "scenario 2 (IHB)",\
     input_case1_kkprime u 1:4:5 with yerrorbars ls 10021 notitle "scenario 1  (LC)",\
     input_case2_kkprime u 1:4:5 with yerrorbars ls 10022 notitle "scenario 2 (IHB)"

unset multiplot
set term wxt
