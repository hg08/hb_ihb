#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

#ST_mean_std_343_MB-pol.txt
simType = 'MB-pol' if len(sys.argv) < 2 else sys.argv[1]
numMolecules = np.array([125, 216, 343, 512, 729, 1000])
meanValues = np.zeros(numMolecules.shape)
stdValues = np.zeros(numMolecules.shape)
stdErrors = np.zeros(numMolecules.shape)
for idx, num in enumerate(numMolecules):
    filename = "output/ST_mean_std_{}_{}.txt".format(num, simType)
    data = np.loadtxt(filename)
    meanValues[idx], stdValues[idx], stdErrors[idx] = data[0], data[1], data[2]

plt.figure(figsize=(6, 2.5))
plt.axhline(y=71.73, linestyle='--', label='71.73 mN/m')
plt.errorbar(numMolecules, meanValues, yerr=stdErrors, fmt='o', color='black', capsize=5, ecolor='black', label=simType)
plt.xlabel('Number of water molecules')
plt.ylabel('Surface tension (mN/m)')
plt.xticks(numMolecules)
plt.tick_params(direction="in", axis='both', top=True)
plt.ylim(40, 100)
plt.legend(frameon=False, ncol=2, loc='upper left')
plt.tight_layout()
plt.savefig('ST_{}.png'.format(simType))
plt.savefig('ST_{}.pdf'.format(simType))
plt.show()
