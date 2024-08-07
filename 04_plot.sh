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
rm -f  tmp
keyword="Final set of parameters"
awk -v var="$keyword" '$0~var{x=NR+5; y=NR+2}(NR<=x && NR >=y){print $1, $3, $4, $5}' double_fit_c2_${system}.log > original_para_double_fit_c2_${system}.dat 
# sort in each set of parameters ("DSC") 
awk -v n=4 '1; NR % n == 0 {print ""; print ""}' original_para_double_fit_c2_${system}.dat > tmp 
awk -v OFS='\t' ' NF<pNF || NR==1 { blockNr++ } { print blockNr, NF, NR, (NF>1 ? $2 : NR), $0; pNF=NF }' tmp | sort -n -k1,1 -k2,2 -k4,4 -k3,3 | cut -f5-  > para_double_fit_c2_${system}.dat 
# Now we can remove the blank lines
sed -i '/^$/d' para_double_fit_c2_${system}.dat

mv para_double_fit_c2_${system}.dat ../2_orientation/output
#prapare the scripts and plot
cp plot_Fig5_template.gp plot_Fig5_${system}.gp
sed -i "s/SYSTEM/${system}/g" plot_Fig5_${system}.gp
gnuplot plot_Fig5_${system}.gp

rm tmp double_fit_c2_${system}.log single_fit_c2_${system}.log 

cd ..
