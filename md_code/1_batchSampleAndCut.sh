for numWater in 125 216 343 512 729 1000
do
	echo "Processing $numWater"
	# Get current path
	inputFile="MB-pol/{}h2o/dump.xyz".format(numWater)
	outputFile="{}h2o-240ps-mbx-pol.xyz".format(numWater)
	./sampleAndCut.py $inputFile $outputFile
done
