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

# --- Step m1: Resample trajectory 
cd m1_resample # (
mkdir -p output
rm -f $trajFile
ln -s ../m2_traj/${trajFile} .
for start_frame in 0 10000 20000 30000 40000 50000
do
subTrajFile=${system}_s${start_frame}.xyz
python resample.py $trajFile ${start_frame} output/$subTrajFile
done
cd .. # ) 
# --- End Step m1

# --- Step 0: Generate surface trajectory
cd 0_prepare # (
for start_frame in 0 10000 20000 30000 40000 50000
do
	# a. Prepare input trajectory
	subTraj=${system}_s${start_frame}
	subTrajFile=${subTraj}.xyz
	rm -f $subTrajFile 
	ln -s ../m1_resample/output/$subTrajFile .

	# b. Prepare input metadata
	inputFile=input_${subTraj}
	cp input_template $inputFile
	sed -i "s/SUBTRAJ/${subTraj}/g" $inputFile
	sed -i "s/SIZEX/${sizeX}/g" $inputFile
	sed -i "s/SIZEY/${sizeY}/g" $inputFile
	sed -i "s/SIZEZ/${sizeZ}/g" $inputFile 

	# c. Run Chandler
	python 1_chandler_fast.py < $inputFile  

	# d. Clean up
	mkdir -p output
	mv *.cube output/
	mv surf_${subTraj}.dat output/  
done
cd .. # )
# --- End Step 0
