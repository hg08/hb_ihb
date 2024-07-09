mkdir -p output
rm -f *png

echo "Processing statistics..."

filePath=output

for case in 1_case1 2_case2
do 
    files=("output/128w-pos-1_s0_${case}_layers_c_at_ref.dat"  
    "output/128w-pos-1_s10000_${case}_layers_c_at_ref.dat"
    "output/128w-pos-1_s20000_${case}_layers_c_at_ref.dat" 
    "output/128w-pos-1_s30000_${case}_layers_c_at_ref.dat" 
    "output/128w-pos-1_s40000_${case}_layers_c_at_ref.dat" 
    "output/128w-pos-1_s50000_${case}_layers_c_at_ref.dat")

    mean_c_at_ref=output/${case}_layers_c_at_ref.dat
    rm -rf $mean_c_at_ref
    awk '
    # For each line in each file
    { 
	    count[FNR]++
	    sum[FNR] += $2
	    sumsq[FNR] += ($2*$2)
    }
    # After processing all files
    END {
        for (i = 1; i <= FNR; i++){
	mean = sum[i] / count[i]
	stddev = sqrt(sumsq[i]/count[i] - (mean * mean))
	printf ("%s %f %f\n", $1, mean, stddev)

	}
    } '  "${files[@]}" > $mean_c_at_ref

    sed -i "s/0.000000/        /g" $mean_c_at_ref
done
