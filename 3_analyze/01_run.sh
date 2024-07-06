mkdir -p output
rm *png
for scenario in 1_case1 2_case2;
do
	echo "Processing ${scenario}..."
	# Create symbolic links to the output files
	for i in `seq 1 1 6`; 
	do 
		# The symbol * is used to match 1.dat and 1.0.dat 
		ln -s ../${scenario}/output/128w_itp_wat_pair_hbacf_h_ihb_${i}*.dat c_${i}.dat; 
	done

	# Extract c(t_f) from files `c_${i}.dat`, where t_f = 0, 1, 2,... ps.
	# ns = 1, dt = 0.04 ps, unit = 1 ps, rows = unit/(dt*ns)
	rowNum=25
	for i in `seq 1 1 6`; 
	do 
		awk -v rowNum=$rowNum "NR%rowNum==1{print}" c_${i}.dat > output/c_${i}_at_ref.dat; 
	done

	# Combine values of special ref. time (1, 2, 5, 10) for all layers.
	cd output
	rm -rf layers_c_at_ref_temp.dat
	touch layers_c_at_ref_temp.dat 
	for i in `seq 1 1 6`;
	do  
		awk 'NR==2 || NR==3 || NR==6 || NR==11{print}' c_${i}_at_ref.dat > ${i}_at_ref.dat
		cat layers_c_at_ref_temp.dat ${i}_at_ref.dat > tmp
		mv tmp layers_c_at_ref_temp.dat
		rm ${i}_at_ref.dat
	done
	rm -rf c_*_at_ref.dat
	cd ..

	rm -rf c_*.dat
	# Format the output file
	./formatting.py output/layers_c_at_ref_temp.dat  output/${scenario}_layers_c_at_ref.dat
	rm output/layers_c_at_ref_temp.dat


	# Fit reactive flux constants for all layers
	kkprime=output/${scenario}_kkprime.dat
	rm -rf $kkprime
	for i in {1..6}
	do
		filePath=../${scenario}/output
		./fit_reactive_flux.py ${filePath}/128w_itp_wat_pair_hbacf_h_ihb_${i}*.dat \
			${filePath}/128w_itp_wat_pair_hbacf_n_ihb_${i}*.dat \
			${filePath}/128w_itp_wat_pair_hbacf_k_ihb_${i}*.dat >> $kkprime
	done

        # Add the THICKNESS
	awk '{print NR, $0}' $kkprime > tmp
	mv tmp $kkprime
done

