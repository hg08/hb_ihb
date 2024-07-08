# 1. Calculate C2(t_f).

# Extract C2(t_f) from files `c2_${i}.dat`, where t_f = 0, 1, 2,... ps.
# ns = 1, dt = 0.04 ps, unit = 1 ps, rows = unit/(dt*ns)
rowNum=25
for i in {1..6}
do
    awk -v rowNum=$rowNum "NR%rowNum==1{print}" ../2_orientation/output/128w_itp_c2_ihb_${i}.0.dat  > output/c2_${i}_at_ref.dat
done

# Combine values of special ref. time (1, 2, 5, 10) for all layers.
cd output
rm -rf layers_c2_at_ref_temp.dat
touch layers_c2_at_ref_temp.dat
for i in {1..6}
do
	awk 'NR==2 || NR==3 || NR==6 || NR==11{print}' c2_${i}_at_ref.dat > ${i}_at_ref.dat
	cat layers_c2_at_ref_temp.dat ${i}_at_ref.dat > tmp
	mv tmp layers_c2_at_ref_temp.dat
	rm ${i}_at_ref.dat
done
rm -rf c2_*_at_ref.dat
cd ..

#rm -rf c2_*.dat
# Format the output file
./formatting.py output/layers_c2_at_ref_temp.dat  output/layers_c2_at_ref.dat
#rm output/layers_c2_at_ref_temp.dat

