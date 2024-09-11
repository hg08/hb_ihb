# the input
arr_d=(1 2 3 4 5 6)
arr_nwat=(125 216)
arr_ref=(1 2 5 10)
arr_calc_type=(1_case1_c 2_case2_c)
nsteps=240

time_step=0.04   #IMPORTANT PARAMATER
num_lines_per_ps=25    # $((1/time_step))

for calc_type in ${arr_calc_type[@]}
do
    for d in ${arr_d[@]}
    do 
        for ref in ${arr_ref[@]}
        do 
            # the out
            out=${nsteps}ps-mbx-pol_${calc_type}_at_ref${ref}_${d}A.dat
            rm -f $out output/$out
            
            # extract the data for thickness=d into the file $out
            for nwat in ${arr_nwat[@]}
            do
                awk -v thickness=$d -v ref=${ref} -v nwater=${nwat} -v nlps=${num_lines_per_ps} -v ct=${calc_type} 'NR==1+ref*nlps {print nwater, $0}' output/${nwat}h2o-${nsteps}ps-mbx-pol_${calc_type}_t_${d}.dat >> $out
            done
            mv  $out output
        done
    done
done
