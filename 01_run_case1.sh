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
numAtoms=$(jq -r '.numAtoms' $JSON_FILE)


echo Processing system $system ...
trajFile=$system.xyz 
    
# --- Step 1_case1 
cd 1_case1 # (
# compile and clean up
rm -rf ihb_sulpizi obj
mkdir -p obj
make 
rm -rf modeule_ihb.mod obj
for start_frame in 0 10000 20000 30000 40000 50000
do 
# a. Prepare input trajectory and parameter file
subTraj=${system}_s${start_frame}
subTrajFile=${subTraj}.xyz
surfTrajFile=surf_${subTraj}.dat
rm -rf $subTrajFile $surfTrajFile 
ln -s ../m1_resample/output/$subTrajFile .
ln -s ../0_prepare/output/$surfTrajFile .
for d in {1..6}
do 
	inputFile=input_${subTraj}_${d}
	cp input_template $inputFile
	sed -i "s/SIZEX/${sizeX}/g" $inputFile
	sed -i "s/SIZEY/${sizeY}/g" $inputFile
	sed -i "s/SIZEZ/${sizeZ}/g" $inputFile
	sed -i "s/SYSTEM/${subTraj}/g" $inputFile
	sed -i "s/SUBTRAJ/${subTraj}/g" $inputFile
	sed -i "s/NUMATOM/${numAtom}/g" $inputFile
	sed -i "s/SURFTRAJFILE/${surfTrajFile}/g" $inputFile
	sed -i "s/THICKNESS/${d}/g" $inputFile
done

# b. Run Chandler
for d in {1..6}
do 
	inputFile=input_${subTraj}_${d}
	./ihb_sulpizi < $inputFile 
done

# c. Clean up
mkdir -p output
mv ${subTraj}_* output/
mv recentered*.xyz output/
done
cd .. # )
# --- End step 1_case1
