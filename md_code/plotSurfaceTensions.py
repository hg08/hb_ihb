#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

#ST_mean_std_343_MB-pol.txt
simType = 'MB-pol' if len(sys.argv) < 2 else sys.argv[1]
label = 'TIP4P/2005' if simType == 'TIP4P2005' else simType
numMolecules = np.array([125, 216, 343, 512, 729, 1000])
meanValues = np.zeros(numMolecules.shape)
stdValues = np.zeros(numMolecules.shape)
stdErrors = np.zeros(numMolecules.shape)
for idx, num in enumerate(numMolecules):
    filename = "output/ST_mean_std_{}_{}.txt".format(num, simType)
    data = np.loadtxt(filename)
    meanValues[idx], stdValues[idx], stdErrors[idx] = data[0], data[1], data[2]

cube_roots = np.cbrt(numMolecules)
scale = 5
offset = -0.2
fontSize = 14
plt.rcParams.update({'font.size': fontSize})  
plt.figure(figsize=(1.05*scale, 1*scale))
plt.axhline(y=71.73, linestyle='--', color='grey', label='Experiment')
plt.errorbar(cube_roots, meanValues, yerr=stdErrors, fmt='o', color='black', capsize=5, ecolor='black', label=label)
plt.xticks(cube_roots, labels=numMolecules)  # Label ticks with the actual cubic numbers
plt.xlabel('Number of water molecules')
plt.ylabel('Surface tension $\gamma$ (mN/m)')
plt.tick_params(direction="in", axis='both', which='both', top=True, right=True)
plt.ylim(40, 100)
plt.legend(frameon=False, ncol=2, loc='upper left')
plt.tight_layout()
plt.savefig('ST_{}.png'.format(simType), dpi=300)
plt.savefig('ST_{}.pdf'.format(simType))
#plt.show()
