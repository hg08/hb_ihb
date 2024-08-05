#!/bin/bash

simType="TIP4P2005"
#simType="MB-pol"
for numWater in 125 216 343 512 729 1000
do
    echo "Processing $numWater"
    if [ $simType == "MB-pol" ]; then
    	inputFile="MB-pol/${numWater}h2o/dump.xyz"
    	outputFile="${numWater}h2o-240ps-mbx-pol.xyz"
    else 
    	inputFile="TIP4P2005/${numWater}w/dump.xyz"
    	outputFile="${numWater}w-240ps-tip4p.xyz"
    fi

    ./sampleAndCut.py $inputFile $outputFile
done
