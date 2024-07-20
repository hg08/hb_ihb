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
simTime=$(jq -r '.simTime' $JSON_FILE)
numFrame=$(jq -r '.numFrames' $JSON_FILE)

# --- Step 1_case1 
echo Processing system $system ...

cd 1_case1 # (
# compile and clean up
rm -rf ihb_sulpizi obj
mkdir -p obj
make 
rm -rf modeule_ihb.mod obj

# prepare input position and surface trajectories 
trajFile=$system.xyz 
surfTrajFile=surf_${system}.dat
rm -rf $trajFile $surfTrajFile
ln -s ../m2_traj/$trajFile .
ln -s ../0_prepare/output/$surfTrajFile .

# By setting the number of sub-trajectories and 
# the time length of each sub-trajectory, we can 
# 'divide' the whole trajectory into several sub-trajectories

numSubTraj=6
subTrajTime=35 # in ps
if [ $subTrajTime -gt $simTime ]; then
	echo "Error: subTrajTime is larger than simTime"
	exit 1
fi
offSetInPs=$(echo "($simTime - $subTrajTime) / ($numSubTraj - 1)" | bc -l)
offSet=$(echo "scale=0; $offSetInPs / $dt" | bc -l)
for (( i=0; i<numSubTraj; i++ ))
do # Loop over sub-trajectories
	# We don't need to really split the trajectory file
	# We just need to set the start and end frame for each sub-trajectory
	start_frame=$(echo "$offSet * $i" | bc -l)
	echo "Sub-trajectory $i start frame: $start_frame"
	end_frame=$(echo "scale=0; 1 + $start_frame + $subTrajTime / $dt" | bc -l)
	echo "Sub-trajectory $i end frame: $end_frame"

	# a. Prepare parameter file
	for d in {1..6}
	do 
		inputFile=input_${system}_s${i}_${d}
		cp input_template $inputFile
		sed -i "s/SIZEX/${sizeX}/g" $inputFile
		sed -i "s/SIZEY/${sizeY}/g" $inputFile
		sed -i "s/SIZEZ/${sizeZ}/g" $inputFile
		sed -i "s/DTINPS/${dt}/g" $inputFile
		sed -i "s/SYSTEM/${system}_s${i}/g" $inputFile
		sed -i "s/WHOLETRAJ/${system}/g" $inputFile
		sed -i "s/FRAMESTART/${start_frame}/g" $inputFile
		sed -i "s/FRAMEEND/${end_frame}/g" $inputFile
		sed -i "s/NUMATOM/${numAtom}/g" $inputFile
		sed -i "s/SURFTRAJFILE/${surfTrajFile}/g" $inputFile
		sed -i "s/THICKNESS/${d}/g" $inputFile
	done
	
	# b. Run Chandler
	for d in {1..6}
	do 
		inputFile=input_${system}_s${i}_${d}
		./ihb_sulpizi < $inputFile 
	done
	
	# c. Clean up
	mkdir -p output
	mv ${system}_* output/
	mv recentered*.xyz output/
done
cd .. # )
# --- End step 1_case1
