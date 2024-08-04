#!/bin/bash

for numWater in 125 216 343 512 729 1000
do
    echo "Processing $numWater"
    # Get current path
    inputFile="MB-pol/${numWater}h2o/dump.xyz"
    outputFile="${numWater}h2o-240ps-mbx-pol.xyz"
    ./sampleAndCut.py $inputFile $outputFile
done
