# Recenter 
mkdir -p output

outfile=sorted_density_OH.dat
awk -f process_and_sort.awk output/density_OH.dat > $outfile

mv $outfile output

