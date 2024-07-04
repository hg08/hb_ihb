#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

# Calculate the coefficients k and k_prime in k=kc +k_prime n.

# Get the command line arguments using sys.argv
if len(sys.argv) != 4:
    print("Usage: python fit_reactive_flux.py <c> <n> <k>")
    sys.exit(1)

c_file = sys.argv[1]
n_file = sys.argv[2]
k_file = sys.argv[3]

#For times beyond the transient period
#IF STUDY LONG-TIME RELAXATION
start = 5
end = 300 

# IF STUDY SHORT-TIME RELAXTION
df_c = pd.read_csv(c_file, header=None)
df_c.columns = ["time","c"]
df_c.head()
#print("shape of c:{}".format(df_c.shape))

# load x2= n(t)
df_n = pd.read_csv(n_file, header=None)
df_n.columns = ["time","n"]
df_n.head()
#print("shape of n:{}".format(df_n.shape))

df_k = pd.read_csv(k_file, header=None)
df_k.columns = ["time","k"]
df_k.head()

c_index= df_c['time'].values
c=df_c['c'].values
c = c[start:end]
c_index = c_index[start:end]

n=df_n['n'].values
n = n[start:end]
#print("n.shape={}".format(n.shape))

k=df_k['k'].values
k = k[start:end]
#print("shape of k:{}".format(df_k.shape))

A = np.stack((c, n), 1)

# If we write $c(t)$, $n(t)$ and $k(t)$ as vectors,and denote 
# them as ${\bf c}$,${\bf n}$ and ${\bf k}$. The $ k$ and $k'$ 
# can be determined from the matrix $ A = ({\bf c} \ {\bf n})$, 
# which is of multiple rows and 2 columns. 
# We can obtain $k$ and $k'$ by
# $[k ,\ -k']^T= (A^T A)^{-1} A^T {\bf k}$.
A_trans = A.T 
B = np.matmul(A_trans, A)
# Calculate the inverse of Matrix B
B_inv = np.linalg.inv(B)
b1 = np.matmul(A_trans,k)
x = np.matmul(B_inv, b1)

plt.plot(c_index,x[0]*c +x[1]*n, '--')
plt.plot(c_index,k, 'r-')
plt.savefig('least_square_fit_reactive_flux_ADH.png')

#print("(k, k') = {} {}".format(x[0], -x[1]))
print("{} {}".format(x[0], -x[1]))
