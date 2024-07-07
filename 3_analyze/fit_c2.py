#!/usr/bin/env python3
#Purpose:
#    Find the decay rate tau2. c2(t) = exp(-t/b).
#    This relation is transformed to a linear fit: lnc2 = -1/b

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys, datetime

# Get the command line arguments using sys.argv
if len(sys.argv) != 2:
    print("Usage: python fit_c2.py <c2>")
    sys.exit(1)

c2_file = sys.argv[1]

#For times [0: end] 
#IF STUDY THE DECAY OF C2.
end = 150 
dt = 0.04

data_c2 = np.loadtxt(c2_file)

simTime = data_c2[:,0] 
timeMask = (simTime < end*dt)
simTime = simTime[timeMask]
c2 = data_c2[timeMask,1]
lnc2 = np.log(c2)

# We write $lnc2(t)$, $t$  as vectors, and denote 
# them as lnc2 and t. The slope = -1/tau2 
# can be determined from the matrix A = (t), 
# which is of multiple rows and 1 columns. 
# We can obtain the slope =-1/tau2 by
# slope = (A^T A)^{-1} A^T lnc2.
A = simTime 
A_trans = A.T
# B is a scalar
B = np.matmul(A_trans, A)
B_inv = 1.0/B
b1 = np.matmul(A_trans, lnc2)
# x is a scalar
x = B_inv * b1
tau2 = -1.0/x
plt.plot(simTime, x*simTime, '--')
plt.plot(simTime, lnc2, 'r-')
# 
#plt.savefig('c2_decay_rate_fit_' + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-%f") + '.png')

print("{}".format(tau2))
