# the input
arr_d=(1 2 3 4 5 6)
arr_nwat=(125 216 343 512)
nsteps=240


for d in ${arr_d[@]}
do 
    # the out
    out=${nsteps}ps-mbx-pol_kkprime_${d}A.dat
    rm -f $out output/$out
    
    # extract the data for thickness=d into the file $out
    for nwat in ${arr_nwat[@]}
    do
        awk -v thickness=$d -v nwater=${nwat} 'NR==thickness {print nwater, $0}' output/${nwat}h2o-${nsteps}ps-mbx-pol_kkprime.dat >> $out
    done
    mv  $out output
done
