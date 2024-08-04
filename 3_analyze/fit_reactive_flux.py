#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys, datetime

# Get the command line arguments using sys.argv
if len(sys.argv) != 4:
    print("Usage: python fit_reactive_flux.py <c> <n> <k>")
    sys.exit(1)

c_file = sys.argv[1]
n_file = sys.argv[2]
k_file = sys.argv[3]

#For times beyond the transient period
#IF STUDY LONG-TIME RELAXATION
window_size = 10  # ps
dt = 0.04 # ps
start = 4 # Get rid of the first few fluctuating data points
end = int(window_size / dt) + 1 

# IF STUDY SHORT-TIME RELAXTION
data_c = np.loadtxt(c_file)
data_n = np.loadtxt(n_file)
data_k = np.loadtxt(k_file)

simTime = data_c[:,0] 
timeMask = (simTime >= start*dt) & (simTime < end*dt)
simTime = simTime[timeMask]
c = data_c[timeMask,1]
n = data_n[timeMask,1]
k = data_k[timeMask,1]

A = np.stack((c, n), 1)
# If we write $c(t)$, $n(t)$ and $k(t)$ as vectors,and denote 
# them as ${\bf c}$,${\bf n}$ and ${\bf k}$. The $ k$ and $k'$ 
# can be determined from the matrix $ A = ({\bf c} \ {\bf n})$, 
# which is of multiple rows and 2 columns. 
# We can obtain $k$ and $k'$ by
# $[k ,\ -k']^T= (A^T A)^{-1} A^T {\bf k}$.
A_trans = A.T 
B = np.matmul(A_trans, A)
B_inv = np.linalg.inv(B)
b1 = np.matmul(A_trans,k)
x = np.matmul(B_inv, b1)
plt.plot(simTime, x[0]*c +x[1]*n, '--')
plt.plot(simTime, k, 'r-')
# 
#plt.savefig('kkprime_fit_' + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-%f") + '.png')

print("{} {}".format(x[0], -x[1]))
