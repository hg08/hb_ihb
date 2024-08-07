# Find the latest json file in m2_traj
if [ $# -eq 0 ]; then
    JSON_FILE=$(ls -t m2_traj/*.json | head -1)
else
    JSON_FILE=$1
fi

system=$(jq -r '.system' $JSON_FILE)
sizeX=$(jq -r '.sizeX' $JSON_FILE)
sizeY=$(jq -r '.sizeY' $JSON_FILE)
sizeZ=$(jq -r '.sizeZ' $JSON_FILE)
numAtom=$(jq -r '.numAtoms' $JSON_FILE)



# 4_plot
cd 4_plot
gnuplot -e "system='$system'" ./plot_Fig2.gp
gnuplot ./plot_Fig3.gp
gnuplot ./plot_Fig4.gp
gnuplot ./plot_Fig5.gp
gnuplot ./plot_FigS2_${system}.gp
gnuplot ./plot_FigS3_${system}.gp

rm -f fit.log
gnuplot ./plot_S5_${system}.gp &&
mv fit.log  single_fit_c2_${system}.log

rm -f fit.log
gnuplot ./plot_S6_${system}.gp &&
mv fit.log  double_fit_c2_${system}.log
# extract the fittting parameters
rm -f tmp0 tmp
keyword="Final set of parameters"
awk -v var="$keyword" '$0~var{x=NR+5; y=NR+2}(NR<=x && NR >=y){print  $1, $3, $4, $5}' double_fit_c2_${system}.log > tmp0 
sort -t, -nk1 tmp0 > tmp 
# The "1;" is a condition that is always true, and will trigger the default action which is to print the current line.
awk -v n=4 '1; NR % n == 0 {print ""; print ""}' tmp > para_double_fit_c2_${system}.dat 

cd ..
