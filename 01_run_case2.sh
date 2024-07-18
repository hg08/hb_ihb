# Find the latest json file in m2_traj
if [ $# -eq 0 ]; then
    JSON_FILE=$(ls -t m2_traj/*.json | head -1)
else
    JSON_FILE=$1
fi

name=$(jq -r '.name' $JSON_FILE)
system=$(jq -r '.system' $JSON_FILE)
dt=$(jq -r '.dt' $JSON_FILE)
sizeX=$(jq -r '.sizeX' $JSON_FILE)
sizeY=$(jq -r '.sizeY' $JSON_FILE)
sizeZ=$(jq -r '.sizeZ' $JSON_FILE)
numAtom=$(jq -r '.numAtoms' $JSON_FILE)
numFrame=$(jq -r '.numFrames' $JSON_FILE)
frameStart=0
frameEnd=$(($numFrame-1))

# --- Step 2_case2
echo Processing system $system ...
trajFile=${name} 
cd 2_case2 # (
# Compile and clean up
./compile.sh
# a. Prepare input trajectory and parameter file
surfTrajFile=surf_${system}.dat
rm -rf $trajFile $surfTrajFile 
ln -s ../m2_traj/$trajFile .
ln -s ../0_prepare/output/$surfTrajFile .
for d in {1..6}
do 
	inputFile=input_${system}_${d}
	cp input_template $inputFile
	sed -i "s/SIZEX/${sizeX}/g" $inputFile
	sed -i "s/SIZEY/${sizeY}/g" $inputFile
	sed -i "s/SIZEZ/${sizeZ}/g" $inputFile
	sed -i "s/DTINPS/${dt}/g" $inputFile
	sed -i "s/SYSTEM/${system}/g" $inputFile
	sed -i "s/FRAMESTART/${frameStart}/g" $inputFile
	sed -i "s/FRAMEEND/${frameEnd}/g" $inputFile
	sed -i "s/NUMATOM/${numAtom}/g" $inputFile
	sed -i "s/SURFTRAJFILE/${surfTrajFile}/g" $inputFile
	sed -i "s/THICKNESS/${d}/g" $inputFile
done

# b. Run ihb 
for d in {1..6}
do 
	inputFile=input_${system}_${d}
	./main_interface_c_format2 < $inputFile
	./main_interface_k_format2 < $inputFile
	./main_interface_n_format2 < $inputFile
done

# c. Clean up
mkdir -p output
mv ${system}_* output/

cd .. # ) 
# --- End step 2_case2
