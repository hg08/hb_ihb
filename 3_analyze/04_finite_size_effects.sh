d=4
arr_nwat=(125 216 343 512)
nsteps=240
rm ${nsteps}ps-mbx-pol_kkprime_${d}A.dat

for nwat in ${arr_nwat[@]}
do
    awk -v thickness=$d -v nwater=${nwat} 'NR==thickness {print nwater, $0}' output/${nwat}h2o-${nsteps}ps-mbx-pol_kkprime.dat >> ${nsteps}ps-mbx-pol_kkprime_${d}A.dat
done
