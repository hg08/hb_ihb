#!/usr/bin/env python
import numpy as np

# Obtain command line arguments
import sys
if len(sys.argv) != 3:
    print("Usage: python formatting.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

data = np.loadtxt(input_file)

# Function to write the data block to a file
def write2file(block, t, f):
    t = int(t)
    f.write("#t={} ps\n".format(t))
    for i in range(block.shape[0]):
        d = i + 1 # This is the layer number
        f.write("{}        {}\n".format(d, block[i, 1])) 
    f.write("\n\n")

# Create a file to write the data
f = open(output_file, 'w')

timeSteps = np.array([1.0, 2.0, 5.0, 10.0])
for t in timeSteps:
    print("Time step: ", t)
    mask = data[:, 0] == t
    data_t = data[mask]
    write2file(data_t, t, f)
f.close()

