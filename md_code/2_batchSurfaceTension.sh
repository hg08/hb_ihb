#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    echo "Usage: 2_batchSurfaceTension.sh <simType>"
    exit 1
fi

simType=$1
#simType="TIP4P2005"
#simType="MB-pol" 
for numWater in 125 216 343 512 729 1000
do
    echo "Processing $simType $numWater"
    if [ $simType == "MB-pol" ]; then
    	inputFile="MB-pol/${numWater}h2o/log.lammps"
    else 
    	inputFile="TIP4P2005/${numWater}w/log.lammps"
    fi

    python surfaceTension.py $inputFile # Output to output folder
done
