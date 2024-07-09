# Systems
systems=("128w-pos-1")
sizeXs=("15.6404")
sizeYs=("15.6404")
sizeZs=("31.2808")
numAtoms=("384")

# Calculation Workflow 
for i in "${!systems[@]}" 
do
    system=${systems[$i]}
    sizeX=${sizeXs[$i]}
    sizeY=${sizeYs[$i]}
    sizeZ=${sizeZs[$i]}
    numAtom=${numAtoms[$i]}

    echo Processing system $system ...
    trajFile=$system.xyz 

	# --- Step m1: Resample trajectory 
	cd m1_resample # (
	mkdir -p output
	rm -f $trajFile
	ln -s ../m2_traj/${trajFile} .
	for start_frame in 0 10000 20000 30000 40000 50000
	do
	subTrajFile=${system}_s${start_frame}.xyz
	python resample.py $trajFile ${start_frame} output/$subTrajFile
	done
	cd .. # ) 
	# --- End Step m1

	# --- Step 0: Generate surface trajectory
	cd 0_prepare # (
	for start_frame in 0 10000 20000 30000 40000 50000
	do
		# a. Prepare input trajectory
		subTraj=${system}_s${start_frame}
		subTrajFile=${subTraj}.xyz
		rm -f $subTrajFile 
		ln -s ../m1_resample/output/$subTrajFile .

		# b. Prepare input metadata
		inputFile=input_${subTraj}
		cp input_template $inputFile
		sed -i "s/SUBTRAJ/${subTraj}/g" $inputFile
		sed -i "s/SIZEX/${sizeX}/g" $inputFile
		sed -i "s/SIZEY/${sizeY}/g" $inputFile
		sed -i "s/SIZEZ/${sizeZ}/g" $inputFile 

		# c. Run Chandler
		python 1_chandler_fast.py < $inputFile  

		# d. Clean up
		mkdir -p output
		mv *.cube output/
		mv surf_${subTraj}.dat output/  
	done
	cd .. # )
	# --- End Step 0
    
	# --- 1_case1 
	cd 1_case1 # (
	# compile and clean up
	rm -rf ihb_sulpizi 
	mkdir -p obj
	make 
	rm -rf modeule_ihb.mod obj
	for start_frame in 0 10000 20000 30000 40000 50000
	do 
	# a. Prepare input trajectory and parameter file
	subTraj=${system}_s${start_frame}
	subTrajFile=${subTraj}.xyz
	surfTrajFile=surf_${subTraj}.dat
	rm -rf $subTrajFile $surfTrajFile 
	ln -s ../m1_resample/output/$subTrajFile .
	ln -s ../0_prepare/output/$surfTrajFile .
	for d in {1..6}
	do 
		inputFile=input_${subTraj}_${d}
		cp input_template $inputFile
		sed -i "s/SIZEX/${sizeX}/g" $inputFile
		sed -i "s/SIZEY/${sizeY}/g" $inputFile
		sed -i "s/SIZEZ/${sizeZ}/g" $inputFile
		sed -i "s/SYSTEM/${subTraj}/g" $inputFile
		sed -i "s/SUBTRAJ/${subTraj}/g" $inputFile
		sed -i "s/NUMATOM/${numAtom}/g" $inputFile
		sed -i "s/SURFTRAJFILE/${surfTrajFile}/g" $inputFile
		sed -i "s/THICKNESS/${d}/g" $inputFile
	done

	# b. Run Chandler
	for d in {1..6}
	do 
		inputFile=input_${subTraj}_${d}
		./ihb_sulpizi < $inputFile 
	done

	# c. Clean up
	mkdir -p output
	mv ${subTraj}_* output/
	mv recentered*.xyz output/
	done
	cd .. # )
	# --- End 1_case1

	# --- 2_case2
	cd 2_case2 # (
	# Compile and clean up
	./compile.sh
	for start_frame in 0 10000 20000 30000 40000 50000
	do 
		# a. Prepare input trajectory and parameter file
		subTraj=${system}_s${start_frame}
		subTrajFile=${subTraj}.xyz
		surfTrajFile=surf_${subTraj}.dat
		rm -rf $subTrajFile $surfTrajFile 
		ln -s ../m1_resample/output/$subTrajFile .
		ln -s ../0_prepare/output/$surfTrajFile .
		for d in {1..6}
		do 
			inputFile=input_${subTraj}_${d}
				cp input_template $inputFile
			sed -i "s/SIZEX/${sizeX}/g" $inputFile
			sed -i "s/SIZEY/${sizeY}/g" $inputFile
			sed -i "s/SIZEZ/${sizeZ}/g" $inputFile
			sed -i "s/SYSTEM/${subTraj}/g" $inputFile
			sed -i "s/SUBTRAJ/${subTraj}/g" $inputFile
			sed -i "s/NUMATOM/${numAtom}/g" $inputFile
			sed -i "s/SURFTRAJFILE/${surfTrajFile}/g" $inputFile
			sed -i "s/THICKNESS/${d}/g" $inputFile
		done

		# b. Run ihb 
		for d in {1..6}
		do 
			inputFile=input_${subTraj}_${d}
			./main_interface_c_format2 < $inputFile
			./main_interface_k_format2 < $inputFile
				./main_interface_n_format2 < $inputFile
		done

		# c. Clean up
		mkdir -p output
		mv ${subTraj}_* output/

	done
	cd .. # ) 
	# --- End 2_case2

    # --- 2_orientation
    cd 2_orientation # (
    # --- Compile and clean up
    rm -f iad_sulpizi
	mkdir -p obj
	make
	rm -rf obj
	rm module_ihb.mod
    for start_frame in 0 10000 20000 30000 40000 50000
    do 
		# a. Prepare input trajectory and parameter file
		subTraj=${system}_s${start_frame}
		subTrajFile=${subTraj}.xyz
		surfTrajFile=surf_${subTraj}.dat
		rm -rf $subTrajFile $surfTrajFile 
		ln -s ../m1_resample/output/$subTrajFile .
		ln -s ../0_prepare/output/$surfTrajFile .
		for d in {1..6}
		do 
			inputFile=input_${subTraj}_${d}
			cp input_template $inputFile
			sed -i "s/SIZEX/${sizeX}/g" $inputFile
			sed -i "s/SIZEY/${sizeY}/g" $inputFile
			sed -i "s/SIZEZ/${sizeZ}/g" $inputFile
			sed -i "s/SYSTEM/${subTraj}/g" $inputFile
			sed -i "s/SUBTRAJ/${subTraj}/g" $inputFile
			sed -i "s/NUMATOM/${numAtom}/g" $inputFile
			sed -i "s/SURFTRAJFILE/${surfTrajFile}/g" $inputFile
			sed -i "s/THICKNESS/${d}/g" $inputFile
		done

		# b. Run 
		for d in {1..6}
		do 
			inputFile=input_${subTraj}_${d}
			./iad_sulpizi < $inputFile
		done
		
		# c. Clean up
		mkdir -p output
		mv ${subTraj}_c2*dat output/
		rm gmon.out
		rm recentered_traj_pos_sampled.xyz
		rm ${subTraj}_*_list.dat
    done
    cd .. # ) 
       # --- End 2_case2
       
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
				ln -s ../${scenario}/output/${subTraj}_wat_pair_hbacf_h_ihb_${d}*.dat c_${d}.dat; 
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
				./fit_reactive_flux.py ${filePath}/${subTraj}_wat_pair_hbacf_h_ihb_${d}*.dat \
					${filePath}/${subTraj}_wat_pair_hbacf_n_ihb_${d}*.dat \
					${filePath}/${subTraj}_wat_pair_hbacf_k_ihb_${d}*.dat >> $kkprime
			done

			# Add the THICKNESS
			awk '{print NR, $0}' $kkprime > tmp
			mv tmp $kkprime
		done

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
		done	


		# c. Fit tau 2 
		echo "(Subtrj: ${subTraj}) Processing C2 decay rate..."

		# Fit decay constant tau2 of C2 for all layers
		tau2=output/${subTraj}_tau2.dat
		rm -rf $tau2
		for d in {1..6}
		do
				filePath=../2_orientation/output
				./fit_c2.py ${filePath}/${subTraj}_c2_ihb_${d}*.dat >> $tau2
		done

		# Add the THICKNESS
		awk '{print NR, $0}' $tau2 > tmp
		mv tmp $tau2


    done
	cd .. # )

	# 3_statistics-c
        mkdir -p 3_analyze/output
        echo "Processing statistics..."
        
	# Statistics for 1_case1
	files1=() # Initialize an empty array, for saving "output/${subTraj}_${scenario}_layers_c_at_ref.dat"
	for start_frame in `seq 0 10000 50000`
	do 
	    subTraj=${system}_s${start_frame}
            files1+=("3_analyze/output/${subTraj}_1_case1_layers_c_at_ref.dat")
        done
        mean_c_at_ref=3_analyze/output/1_case1_layers_c_at_ref.dat
        rm -rf $mean_c_at_ref
        awk '
        # For each line in each file
	{
		if (FNR > maxFNR) {
		    maxFNR = FNR;
		}
		ndex[FNR] = $1;
                count[FNR]++
                sum[FNR] += $2
                sumsq[FNR] += ($2*$2)
        }

        # After processing all files
        END {
            for (i = 1; i <= FNR; i++){
            mean = sum[i] / count[i]
            stddev = sqrt(sumsq[i]/count[i] - (mean * mean))
	    std_err = stddev/ sqrt(count[i])
            printf ("%s %f %f\n", ndex[i], mean, std_err)
        
            }
        } '  "${files1[@]}" > $mean_c_at_ref
        
        sed -i "s/0.000000/        /g" $mean_c_at_ref
	# END--Statistics for 1_case1
	# Statistics for 2_case2
	files2=() # Initialize an empty array, for saving "output/${subTraj}_${scenario}_layers_c_at_ref.dat"
	for start_frame in `seq 0 10000 50000`
	do 
	    subTraj=${system}_s${start_frame}
            files2+=("3_analyze/output/${subTraj}_2_case2_layers_c_at_ref.dat")
        done
        mean_c_at_ref=3_analyze/output/2_case2_layers_c_at_ref.dat
        rm -rf $mean_c_at_ref
        awk '
        # For each line in each file
	{
		if (FNR > maxFNR) {
		    maxFNR = FNR;
		}
		ndex[FNR] = $1;
                count[FNR]++
                sum[FNR] += $2
                sumsq[FNR] += ($2*$2)
        }
        # After processing all files
        END {
            for (i = 1; i <= FNR; i++){
            mean = sum[i] / count[i]
            stddev = sqrt(sumsq[i]/count[i] - (mean * mean))
	    std_err = stddev/ sqrt(count[i])
            printf ("%s %f %f\n", ndex[i], mean, std_err)
        
            }
        } '  "${files2[@]}" > $mean_c_at_ref
        
        sed -i "s/0.000000/        /g" $mean_c_at_ref
	# END--Statistics for 2_case2
	# END---3_statistics-c

	# 3_statistics-kkprime
        mkdir -p 3_analyze/output
        echo "Processing statistics kkprime"
        
	# Statistics for 1_case1 kkprime
	files1=() # Initialize an empty array, for saving "output/${subTraj}_1_case1_kkprime.dat"
	for start_frame in `seq 0 10000 50000`
	do 
	    subTraj=${system}_s${start_frame}
            files1+=("3_analyze/output/${subTraj}_1_case1_kkprime.dat")
        done
        mean_kkprime=3_analyze/output/1_case1_kkprime.dat
        rm -rf $mean_kkprime
        awk '
        # For each line in each file
	{
		if (FNR > maxFNR) {
		    maxFNR = FNR;
		}
		ndex[FNR] = $1;
                count[FNR]++
                sum[FNR] += $2
                sumsq[FNR] += ($2*$2)
                sum_prime[FNR] += $3
                sumsq_prime[FNR] += ($3*$3)
        }

        # After processing all files
        END {
            for (i = 1; i <= FNR; i++){
            mean = sum[i] / count[i]
            stddev = sqrt(sumsq[i]/count[i] - (mean * mean))
	    std_err = stddev/ sqrt(count[i])
            mean_prime = sum_prime[i] / count[i]
            stddev_prime = sqrt(sumsq_prime[i]/count[i] - (mean_prime * mean_prime))
	    std_err_prime = stddev_prime/ sqrt(count[i])
            printf ("%s %f %f %f %f\n", ndex[i], mean, std_err, mean_prime, std_err_prime)
        
            }
        } '  "${files1[@]}" > $mean_kkprime
        
	# END--Statistics for 1_case1 kkprime
	# Statistics for 2_case2 kkprime
	files2=() # Initialize an empty array, for saving "output/${subTraj}_2_case2_kkprime.dat"
	for start_frame in `seq 0 10000 50000`
	do 
	    subTraj=${system}_s${start_frame}
            files2+=("3_analyze/output/${subTraj}_2_case2_kkprime.dat")
        done
        mean_kkprime=3_analyze/output/2_case2_kkprime.dat
        rm -rf $mean_kkprime
        awk '
        # For each line in each file
	{
		if (FNR > maxFNR) {
		    maxFNR = FNR;
		}
		ndex[FNR] = $1;
                count[FNR]++
                sum[FNR] += $2
                sumsq[FNR] += ($2*$2)
                sum_prime[FNR] += $3
                sumsq_prime[FNR] += ($3*$3)
        }

        # After processing all files
        END {
            for (i = 1; i <= FNR; i++){
            mean = sum[i] / count[i]
            stddev = sqrt(sumsq[i]/count[i] - (mean * mean))
	    std_err = stddev/ sqrt(count[i])
            mean_prime = sum_prime[i] / count[i]
            stddev_prime = sqrt(sumsq_prime[i]/count[i] - (mean_prime * mean_prime))
	    std_err_prime = stddev_prime/ sqrt(count[i])
            printf ("%s %f %f %f %f\n", ndex[i], mean, std_err, mean_prime, std_err_prime)
        
            }
        } '  "${files2[@]}" > $mean_kkprime
	# END--Statistics for 2_case2 kkprime
	# END---3_statistics-kkprime
	# 3_statistics-c2
        mkdir -p 3_analyze/output
        echo "Processing statistics c2..."
        
	files=() # Initialize an empty array, for saving "output/*_layers_c2_at_ref.dat"
	for start_frame in `seq 0 10000 50000`
	do 
	    subTraj=${system}_s${start_frame}
            files+=("3_analyze/output/${subTraj}_layers_c2_at_ref.dat")
        done
        mean_c2=3_analyze/output/layers_c2_at_ref.dat
        rm -rf $mean_c2
        awk '
        # For each line in each file
	{
		if (FNR > maxFNR) {
		    maxFNR = FNR;
		}
		ndex[FNR] = $1;
                count[FNR]++
                sum[FNR] += $2
                sumsq[FNR] += ($2*$2)
        }

        # After processing all files
        END {
            for (i = 1; i <= FNR; i++){
            mean = sum[i] / count[i]
            stddev = sqrt(sumsq[i]/count[i] - (mean * mean))
	    std_err = stddev/ sqrt(count[i])
            printf ("%s %f %f\n", ndex[i], mean, std_err)
        
            }
        } '  "${files[@]}" > $mean_c2
        sed -i "s/0.000000/        /g" $mean_c2
   
	# END--3_statistics-c2
	# 3_statistics-tau2
        echo "Processing statistics tau2..."
	files=() # Initialize an empty array, for saving "output/${subTraj}_tau2.dat"
	for start_frame in `seq 0 10000 50000`
	do 
	    subTraj=${system}_s${start_frame}
            files+=("3_analyze/output/${subTraj}_tau2.dat")
        done
        mean_tau2=3_analyze/output/tau2.dat
        rm -rf $mean_tau2
        awk '
        # For each line in each file
	{
		if (FNR > maxFNR) {
		    maxFNR = FNR;
		}
		ndex[FNR] = $1;
                count[FNR]++
                sum[FNR] += $2
                sumsq[FNR] += ($2*$2)
        }

        # After processing all files
        END {
            for (i = 1; i <= FNR; i++){
            mean = sum[i] / count[i]
            stddev = sqrt(sumsq[i]/count[i] - (mean * mean))
	    std_err = stddev/ sqrt(count[i])
            printf ("%s %f %f\n", ndex[i], mean, std_err)
        
            }
        } '  "${files[@]}" > $mean_tau2
	# END---3_statistics-tau2
        

	# 4_plot
	cd 4_plot
	    gnuplot ./plot_Fig3.gp
	    gnuplot ./plot_Fig4.gp
	    gnuplot ./plot_Fig5.gp
	    gnuplot ./plot_Fig7.gp
	cd ..
done
