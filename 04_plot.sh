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
cd ..
