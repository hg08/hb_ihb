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

       
    # --- 3_statistics-c
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
    
    sed -i "s/0.000000/        /g" $mean_c_at_ref ## WARNING: 0.000000 is not general.
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
    

    # 3_statistics-C(t)
    for scenario in 1_case1 2_case2
    do 
        for d in {1..6} 
        do 
            files=() # Initialize an empty array, for saving "output/${subTraj}_${scenario}_**_h_*.dat"
            for start_frame in `seq 0 10000 50000` 
            do
                subTraj=${system}_s${start_frame}
    	        files+=(${scenario}/output/${subTraj}_wat_pair_hbacf_h_ihb_${d}.dat)
                # The symbol * is used to match 1.dat and 1.0.dat 
            done
            mean_c_t=3_analyze/output/${system}_${scenario}_c_t_${d}.dat

    	    rm -rf $mean_c_t
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
            } '  "${files[@]}" > $mean_c_t

        done
    done
    # END--3_statistics-C(t)



	# -- 3_statistics-kkprime
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

	# --- 3_statistics-c2
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
	    gnuplot ./plot_Fig2.gp
	    gnuplot ./plot_Fig3.gp
	    gnuplot ./plot_Fig4.gp
	    gnuplot ./plot_Fig5.gp
	cd ..
#done
