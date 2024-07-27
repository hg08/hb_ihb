# -------------------------
# This script prepares the
# surface trajectory
# -------------------------

# Find the latest json file in m2_traj
if [ $# -eq 0 ]; then
    JSON_FILE=$(ls -t m2_traj/*.json | head -1)
else
    JSON_FILE=$1
fi

# Obtain system information
name=$(jq -r '.name' $JSON_FILE)
system=$(jq -r '.system' $JSON_FILE)
sizeX=$(jq -r '.sizeX' $JSON_FILE)
sizeY=$(jq -r '.sizeY' $JSON_FILE)
sizeZ=$(jq -r '.sizeZ' $JSON_FILE)
numAtom=$(jq -r '.numAtoms' $JSON_FILE)
numFrame=$(jq -r '.numFrames' $JSON_FILE)


# --- Step 0: Generate surface trajectory
echo Processing system $system ...
cd 0_prepare # (
# a. Prepare input trajectory
trajFile=${name} 
rm -f $trajFile
ln -s ../m2_traj/$trajFile .

# b. Prepare input metadata
inputFile=input_${system}
cp input_template $inputFile
sed -i "s/SUBTRAJ/${system}/g" $inputFile
sed -i "s/SIZEX/${sizeX}/g" $inputFile
sed -i "s/SIZEY/${sizeY}/g" $inputFile
sed -i "s/SIZEZ/${sizeZ}/g" $inputFile 
sed -i "s/NUMFRAME/${numFrame}/g" $inputFile

# c. Run Chandler
#python 1_chandler_fast.py < $inputFile  
python 2_chandler_super_fast.py < $inputFile  

# d. Clean up
mkdir -p output
mv *.cube output/
mv surf_${system}.dat output/  
cd .. # )
# --- End Step 0
