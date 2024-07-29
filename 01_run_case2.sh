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
numSubTrajCase1=$(jq -r '.numSubTrajCase1' $JSON_FILE)
subTrajTimeCase1=$(jq -r '.subTrajTimeCase1' $JSON_FILE)
numSubTrajCase2=$(jq -r '.numSubTrajCase2' $JSON_FILE)
subTrajTimeCase2=$(jq -r '.subTrajTimeCase2' $JSON_FILE)
numSubTrajOri=$(jq -r '.numSubTrajOri' $JSON_FILE)
subTrajTimeOri=$(jq -r '.subTrajTimeOri' $JSON_FILE)

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

numSubTraj=$numSubTrajCase2
subTrajTime=$subTrajTimeCase2 # in ps
if (( $(echo "$subTrajTime > $simTime" | bc -l) )); then
        echo "Error: subTrajTime is larger than simTime"
        exit 1
fi

if [ $numSubTraj -eq 1 ]; then
	offSetInPs=0
else
	offSetInPs=$(echo "($simTime - $subTrajTime) / ($numSubTraj - 1)" | bc -l)
fi

offSet=$(echo "scale=0; $offSetInPs / $dt" | bc -l)
for (( i=0; i<numSubTraj; i++ ))
do # Loop over sub-trajectories
        start_frame=$(echo "$offSet * $i" | bc -l)
        echo "Sub-trajectory $i start frame: $start_frame"
        #end_frame=$(echo "scale=0; 1 + $start_frame + $subTrajTime / $dt" | bc -l)
        end_frame=$(echo "scale=0;  $start_frame + $subTrajTime / $dt" | bc -l)
        echo "Sub-trajectory $i end frame: $end_frame"
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
	
	# b. Run ihb 
	for d in {1..6}
	do 
		inputFile=input_${system}_s${i}_${d}
		./main_interface_c_format2 < $inputFile
		./main_interface_k_format2 < $inputFile
		./main_interface_n_format2 < $inputFile
	done

	# c. Clean up
	mkdir -p output
	mv ${system}_* output/

done # End loop over sub-trajectories
cd .. # ) 
# --- End step 2_case2
