rm -rf output--TEST

# Run the test program
python3 1_chandler_fast.py <  input_128w-pos-1_s0

# Move output files to output--TEST directory
mkdir -p output--TEST
mv *.cube output--TEST
