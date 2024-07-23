for numWater in 125 216 343 512 729 1000
do
	echo "Processing $numWater"
	# Get current path
	basePath=$(pwd)
	path="AirWaterInterface${numWater}/run_TIP4P2005_300_hdf5/run_TIP4P2005_300/"
	cd $path
	reoh dump.xyz # Reorange the order of the atoms
	$basePath/sampleAndCut.py dump_RROH.xyz ${numWater}w-TIP4P2005.xyz
	cd $basePath
done
