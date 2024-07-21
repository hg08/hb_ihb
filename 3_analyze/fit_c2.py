#!/usr/bin/env python3
#Purpose:
#    Find the decay rate tau2. c2(t) = exp(-t/b).
#    This relation is transformed to a linear fit: lnc2 = -1/b

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import sys, datetime

# Get the command line arguments using sys.argv
if len(sys.argv) != 2:
    print("Usage: python fit_c2.py <c2>")
    sys.exit(1)

c2_file = sys.argv[1]


def exp_decay(t,tau):
    return np.exp(-t / tau)


#For times [0: end] 
#IF STUDY THE DECAY OF C2.
end = 150 
dt = 0.04

data_c2 = np.loadtxt(c2_file)

simTime = data_c2[:,0] 
timeMask = (simTime < end*dt)
t_data = simTime[timeMask]
y_data = data_c2[timeMask,1]


# Use curve_fit to fit the data
initial_guess = [1]
params, covariance = curve_fit(exp_decay, t_data, y_data, p0=initial_guess)

# Extract fitted parameter and plot the fit
tau_fitted = params[0]
fitted_y_data = exp_decay(t_data, tau_fitted)

plt.scatter(t_data, y_data, label='Data')
plt.plot(t_data, fitted_y_data, label=f'Fitted function(tau ={tau_fitted:.2f})', color='red')
plt.legend()
plt.xlabel('t')
plt.xlabel('c2')
plt.title('Exponential Decay Fit')

# 
#plt.savefig('c2_decay_rate_fit_' + datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-%f") + '.png')

print("{}".format(tau_fitted))
