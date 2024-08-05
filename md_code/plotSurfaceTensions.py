#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

#ST_mean_std_343_MB-pol.txt
simType = 'MB-pol' if len(sys.argv) < 2 else sys.argv[1]
numMolecules = np.array([125, 216, 343, 512, 729, 1000])
meanValues = np.zeros(numMolecules.shape)
stdValues = np.zeros(numMolecules.shape)
for idx, num in enumerate(numMolecules):
    filename = "output/ST_mean_std_{}_{}.txt".format(num, simType)
    data = np.loadtxt(filename)
    meanValues[idx], stdValues[idx] = data[0], data[1]

plt.errorbar(numMolecules, meanValues, yerr=stdValues, fmt='o', color='black', capsize=5, ecolor='black', elinewidth=2, capthick=2, label=simType)
plt.axhline(y=71.73, color='k', linestyle='--', label='Surface Tension Exp. 300 K: 71.73 mN/m')
plt.xlabel('Number of water molecules')
plt.ylabel('Surface tension (mN/m)')
plt.ylim(30, 110)
plt.legend()
plt.savefig('ST_{}.png'.format(simType))
plt.show()
