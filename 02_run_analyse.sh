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

# # Systems
# systems=("128w-pos-1")
# sizeXs=("15.6404")
# sizeYs=("15.6404")
# sizeZs=("31.2808")
# numAtoms=("384")
# 
# # Calculation Workflow 
# for i in "${!systems[@]}" 
# do
#     system=${systems[$i]}
#     sizeX=${sizeXs[$i]}
#     sizeY=${sizeYs[$i]}
#     sizeZ=${sizeZs[$i]}
#     numAtom=${numAtoms[$i]}

    echo Processing system $system ...
    trajFile=$system.xyz 

       
    # --- 3_analyze
    cd 3_analyze # (
    mkdir -p output
    rm -f *.png
    for start_frame in 0 10000 20000 30000 40000 50000
	do 
		subTraj=${system}_s${start_frame}
		# a. Calculate k, k_prime
		for scenario in 1_case1 2_case2
		do
			echo "Processing ${start_frame} ${scenario}..."
			# Create symbolic links to the output files
			for d in {1..6} 
			do 
				# The symbol * is used to match 1.dat and 1.0.dat 
				ln -s ../${scenario}/output/${subTraj}_wat_pair_hbacf_h_ihb_${d}.dat c_${d}.dat; 
			done

			# Extract c(t_f) from files `c_${i}.dat`, where t_f = 0, 1, 2,... ps.
			# ns = 1, dt = 0.04 ps, unit = 1 ps, rows = unit/(dt*ns)
			rowNum=25
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
		done
	done  # End start_frame loop

	# b. Calculate c2
	for start_frame in 0 10000 20000 30000 40000 50000
	do 
		subTraj=${system}_s${start_frame}
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

		#rm -rf c2_*.dat
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

#done
