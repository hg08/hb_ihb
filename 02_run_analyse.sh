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


echo Processing system $system ...
trajFile=$system.xyz 

      
# --- 3_analyze
# a. Calculate k, k_prime
cd 3_analyze # (
mkdir -p output
rm -f *.png

# Case 1
scenario=1_case1
numSubTraj=6
subTrajTime=20 # in ps
if [ $subTrajTime -gt $simTime ]; then
	echo "Error: subTrajTime is larger than simTime"
	exit 1
fi
offSetInPs=$(echo "($simTime - $subTrajTime) / ($numSubTraj - 1)" | bc -l)
offSet=$(echo "scale=0; $offSetInPs / $dt" | bc -l)
for (( i=0; i<numSubTraj; i++ ))
do 
	#start_frame=$(echo "$offSet * $i" | bc -l)
	#cho "Sub-trajectory $i start frame: $start_frame"
	#nd_frame=$(echo "scale=0; 1 + $start_frame + $subTrajTime / $dt" | bc -l)
	#cho "Sub-trajectory $i end frame: $end_frame"

	subTraj=${system}_s${i}
       	echo "Processing ${scenario} ${subTraj} ..."
       	for d in {1..6} 
       	do 
       		# The symbol * is used to match 1.dat and 1.0.dat 
       		ln -s ../${scenario}/output/${subTraj}_wat_pair_hbacf_h_ihb_${d}.dat c_${d}.dat; 
       	done

       	# Extract c(t_f) from files `c_${i}.dat`, where t_f = 0, 1, 2,... ps.
       	# ns = 1, dt = 0.04 ps, unit = 1 ps, rows = unit/(dt*ns)
       	rowNum=25 # TO BE UPDATED
       	for d in {1..6} 
       	do 
       		awk -v rowNum=$rowNum "NR%rowNum==1{print}" c_${d}.dat > output/c_${d}_at_ref.dat; 
       	done

       	# Combine values of special ref. time (1, 2, 5, 10) for all layers.
       	cd output
       	rm -rf layers_c_at_ref_temp.dat
       	touch layers_c_at_ref_temp.dat 
       	for d in {1..6}
       	do  
       		awk 'NR==2 || NR==3 || NR==6 || NR==11{print}' c_${d}_at_ref.dat > ${d}_at_ref.dat
       		cat layers_c_at_ref_temp.dat ${d}_at_ref.dat > tmp
       		mv tmp layers_c_at_ref_temp.dat
       		rm ${d}_at_ref.dat
       	done
       	rm -rf c_*_at_ref.dat
       	cd ..

       	rm -rf c_*.dat
       	# Format the output file
       	./formatting.py output/layers_c_at_ref_temp.dat  output/${subTraj}_${scenario}_layers_c_at_ref.dat
       	rm output/layers_c_at_ref_temp.dat

               
       	# Fit reactive flux constants for all layers
       	kkprime=output/${subTraj}_${scenario}_kkprime.dat
       	rm -rf $kkprime
       	for d in {1..6}
       	do
       		filePath=../${scenario}/output
       		./fit_reactive_flux.py ${filePath}/${subTraj}_wat_pair_hbacf_h_ihb_${d}.dat \
       			${filePath}/${subTraj}_wat_pair_hbacf_n_ihb_${d}.dat \
       			${filePath}/${subTraj}_wat_pair_hbacf_k_ihb_${d}.dat >> $kkprime
       	done

       	# Add the THICKNESS
       	awk '{print NR, $0}' $kkprime > tmp
       	mv tmp $kkprime
done  # End start_frame loop



# Case 2
scenario=2_case2
subTraj=${system} # Only one 'sub-trajectory' for case 2
echo "Processing ${scenario} ${subTraj} ..."
for d in {1..6} 
do 
	# The symbol * is used to match 1.dat and 1.0.dat 
	ln -s ../${scenario}/output/${subTraj}_wat_pair_hbacf_h_ihb_${d}.dat c_${d}.dat; 
done

# Extract c(t_f) from files `c_${i}.dat`, where t_f = 0, 1, 2,... ps.
# ns = 1, dt = 0.04 ps, unit = 1 ps, rows = unit/(dt*ns)
rowNum=25 # TO BE UPDATED
for d in {1..6} 
do 
	awk -v rowNum=$rowNum "NR%rowNum==1{print}" c_${d}.dat > output/c_${d}_at_ref.dat; 
done

# Combine values of special ref. time (1, 2, 5, 10) for all layers.
cd output
rm -rf layers_c_at_ref_temp.dat
touch layers_c_at_ref_temp.dat 
for d in {1..6}
do  
	awk 'NR==2 || NR==3 || NR==6 || NR==11{print}' c_${d}_at_ref.dat > ${d}_at_ref.dat
	cat layers_c_at_ref_temp.dat ${d}_at_ref.dat > tmp
	mv tmp layers_c_at_ref_temp.dat
	rm ${d}_at_ref.dat
done
rm -rf c_*_at_ref.dat
cd ..

rm -rf c_*.dat
# Format the output file
./formatting.py output/layers_c_at_ref_temp.dat  output/${subTraj}_${scenario}_layers_c_at_ref.dat
rm output/layers_c_at_ref_temp.dat

       
# Fit reactive flux constants for all layers
kkprime=output/${subTraj}_${scenario}_kkprime.dat
rm -rf $kkprime
for d in {1..6}
do
	filePath=../${scenario}/output
	./fit_reactive_flux.py ${filePath}/${subTraj}_wat_pair_hbacf_h_ihb_${d}.dat \
		${filePath}/${subTraj}_wat_pair_hbacf_n_ihb_${d}.dat \
		${filePath}/${subTraj}_wat_pair_hbacf_k_ihb_${d}.dat >> $kkprime
done

# Add the THICKNESS
awk '{print NR, $0}' $kkprime > tmp
mv tmp $kkprime



# b. Orientation: calculate c2
echo "Processing orientation..."
for (( i=0; i<numSubTraj; i++ ))
do 
	subTraj=${system}_s${i}
	# Extract C2(t_f) from files `c2_${i}.dat`, where t_f = 0, 1, 2,... ps.
	# ns = 1, dt = 0.04 ps, unit = 1 ps, rows = unit/(dt*ns)
	rowNum=25
	for d in {1..6}
	do
		awk -v rowNum=$rowNum "NR%rowNum==1{print}" ../2_orientation/output/${subTraj}_c2_ihb_${d}.0.dat  > output/c2_${d}_at_ref.dat
	done

       	# Combine values of special ref. time (1, 2, 5, 10) for all layers.
       	cd output
       	rm -rf layers_c2_at_ref_temp.dat
       	touch layers_c2_at_ref_temp.dat
       	for d in {1..6}
       	do
       		awk 'NR==2 || NR==3 || NR==6 || NR==11{print}' c2_${d}_at_ref.dat > ${d}_at_ref.dat
       		cat layers_c2_at_ref_temp.dat ${d}_at_ref.dat > tmp
       		mv tmp layers_c2_at_ref_temp.dat
       		rm ${d}_at_ref.dat
       	done
       	rm -rf c2_*_at_ref.dat
       	cd ..

       	# Format the output file
       	./formatting.py output/layers_c2_at_ref_temp.dat  output/${subTraj}_layers_c2_at_ref.dat
       	rm output/layers_c2_at_ref_temp.dat	

       	# c. Fit tau 2 
       	echo "(Subtrj: ${subTraj}) Processing C2 decay rate..."

       	# Fit decay constant tau2 of C2 for all layers
       	tau2=output/${subTraj}_tau2.dat
       	rm -rf $tau2
       	for d in {1..6}
       	do
       		filePath=../2_orientation/output
       		./fit_c2.py ${filePath}/${subTraj}_c2_ihb_${d}.0.dat >> $tau2
       	done

       	# Add the THICKNESS
       	awk '{print NR, $0}' $tau2 > tmp
       	mv tmp $tau2
done
cd .. # )
