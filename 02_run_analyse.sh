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


echo Processing system $system ...
trajFile=$system.xyz 

      
# --- 3_analyze
# a. Calculate k, k_prime
cd 3_analyze # (
mkdir -p output
rm -f *.png

# Case 1
scenario=1_case1
numSubTraj=$numSubTrajCase1
#subTrajTime=35 # in ps
#if [ $subTrajTime -gt $simTime ]; then
#	echo "Error: subTrajTime is larger than simTime"
#	exit 1
#fi
#offSetInPs=$(echo "($simTime - $subTrajTime) / ($numSubTraj - 1)" | bc -l)
#offSet=$(echo "scale=0; $offSetInPs / $dt" | bc -l)
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
numSubTraj=$numSubTrajCase2
for (( i=0; i<numSubTraj; i++ ))
do
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


# b. Orientation: calculate c2
echo "Processing orientation..."
numSubTraj=$numSubTrajOri
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

# b2. orientation: for specific d, calculate mean of c2
for d in {1..6}
do
    # b2a: combine all c2 files for specific d into one file
    for (( i=0; i<numSubTraj; i++ ))
    do
	subTraj=${system}_s${i}
        rm -f fc2_${i}
    	ln -s 2_orientation/output/${subTraj}_c2_ihb_${d}.0.dat  fc2_${i}
    done
    # b2b: calculate the mean c2 and std error over different subTraj's.
    # obtain the mean c2
    awk '{a[FNR]+=$2;stderr[FNR]+=$2*$2;b[FNR]++;c[FNR]+=$1}END{for(i=1;i<=FNR;i++)print c[i]/b[i],a[i]/b[i],(stderr[i]/b[i]-(a[i]/b[i])**2)/sqrt(b[i]) }' fc2_* > ${system}_c2_ihb_ave_and_stderr_${d}.dat       
    rm fc2_*
    mv ${system}_c2_ihb_ave_*dat 2_orientation/output/
done
# b2c: save info for plotting and fitting the mean c2, and for plotting others
orientation_info_txt=${system}_orientation_info.txt
touch ${orientation_info_txt}
rm -f ${orientation_info_txt}
for d in {1..6}
do 
    echo "input_plot${d} = '../2_orientation/output/${system}_c2_ihb_ave_and_stderr_${d}.dat'" >> ${orientation_info_txt}
done
echo "output_plot = 'c2_${system}.eps'" >> ${orientation_info_txt}
echo "output_plot_single_fit = 'c2_${system}_single_fit.eps'" >> ${orientation_info_txt}
echo "output_plot_double_fit = 'c2_${system}_double_fit.eps'" >> ${orientation_info_txt}
echo "input_plot_fit_para = '../2_orientation/output/para_double_fit_c2_${system}.dat'" >> ${orientation_info_txt}
echo "output_plot_fit_para = 'para_c2_${system}_double_fit.eps'" >> ${orientation_info_txt}
echo "input_case1_kkprime = '../3_analyze/output/${system}_1_case1_kkprime.dat'" >> ${orientation_info_txt}
echo "input_case2_kkprime = '../3_analyze/output/${system}_2_case2_kkprime.dat'" >> ${orientation_info_txt}
echo "input_kkprime = '../3_analyze/output/${system}_kkprime.dat'" >> ${orientation_info_txt}
echo "output_k = 'k_${system}.eps'" >> ${orientation_info_txt}
echo "output_kprime = 'kprime_${system}.eps'" >> ${orientation_info_txt}
# For cHB
for d in {1..6}
do 
echo "input_case1_cHB${d} = '../3_analyze/output/${system}_1_case1_c_t_${d}.dat'" >> ${orientation_info_txt}
echo "input_case2_cHB${d} = '../3_analyze/output/${system}_2_case2_c_t_${d}.dat'" >> ${orientation_info_txt}
done
echo "output_cHB = 'cHB_${system}.eps'" >> ${orientation_info_txt}
# For cHB_at_ref
echo "input_case1_cHB_ref = '../3_analyze/output/${system}_1_case1_layers_c_at_ref.dat'" >> ${orientation_info_txt}
echo "input_case2_cHB_ref = '../3_analyze/output/${system}_2_case2_layers_c_at_ref.dat'" >> ${orientation_info_txt}
echo "output_cHB_ref = 'cHB_ref_${system}.eps'" >> ${orientation_info_txt}
# For c2_at_ref
echo "input_c2_ref = '../3_analyze/output/${system}_layers_c2_at_ref.dat'" >> ${orientation_info_txt}

mv ${orientation_info_txt} 2_orientation/output
# b2d: fit the mean c2


# START density
# Calculate the density of water interface

cd 2_density

#Traj file
rm $system.xyz 
ln -s ../m2_traj/$system.xyz .
#Compile
rm density
gfortran -o density density.f95

# Run and obtain density
stepI=0
stepF=$((stepI + numFrame))
ndiv=501
input_density=input_${system}_density_OH
cp input_density_OH_template $input_density
sed -i "s/SYSTEM/$system/" $input_density
sed -i "s/STEPI/$stepI/" $input_density
sed -i "s/STEPF/$stepF/" $input_density
sed -i "s/NDIV/$ndiv/" $input_density
sed -i "s/SIZEX/$sizeX/" $input_density
sed -i "s/SIZEY/$sizeY/" $input_density
sed -i "s/SIZEZ/$sizeZ/" $input_density


./density < $input_density
mkdir -p output
out_density=density_OH_${system}.dat # This out_density is defined in the "input_density_OH"
mv $out_density output

density_info_txt=${system}_density_info.txt
touch ${density_info_txt}
echo "input_plot = '../2_density/output/${out_density}'" >> ${density_info_txt}
echo "output_plot = 'density_OH_${system}.eps'" >> ${density_info_txt}
echo "sizeZ = '$sizeZ'" >> ${density_info_txt}

mv ${density_info_txt} output

## Recenter (Only for the 128w AIDM simulation )
#outfile=sorted_density_OH.dat
#awk -f process_and_sort.awk output/$out_density > $outfile
#mv $outfile output

#Go back the main direcotry
cd .. 
# END density
