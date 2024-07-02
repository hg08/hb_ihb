# Prepare position file
rm -rf 128w-pos-1_sampled.xyz
ln -s ../m1_resample/128w-pos-1_sampled.xyz 128w-pos-1_sampled.xyz

# Run the program
python3 1_chandler_fast.py < input_chandler_128w

# Move output files to output directory
mkdir -p output
mv *.cube output/
