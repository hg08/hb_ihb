for numWater in 125 216 343 512 729 1000
do
	echo "Processing $numWater"
	# Get current path
	basePath=$(pwd)
	path="MB-pol/125h2o/"
	cd $path
	$basePath/sampleAndCut.py dump.xyz ${numWater}h2o-240ps-mbx-pol.xyz
	cd $basePath
done
