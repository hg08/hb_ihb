#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Plot the surface tension values
# with standard error bars
# ----------------------------
simType = 'MB-pol' if len(sys.argv) < 2 else sys.argv[1]
label = 'TIP4P/2005' if simType == 'TIP4P2005' else simType
numMolecules = np.array([125, 216, 343, 512, 729, 1000])
meanValues = np.zeros(numMolecules.shape)
stdValues = np.zeros(numMolecules.shape)
stdErrors = np.zeros(numMolecules.shape)
n_samples = np.zeros(numMolecules.shape)
for idx, num in enumerate(numMolecules):
    filename = "output/ST_mean_std_{}_{}.txt".format(num, simType)
    data = np.loadtxt(filename)
    meanValues[idx], stdValues[idx], stdErrors[idx], n_samples[idx] = data[0], data[1], data[2], data[3]

cube_roots = np.cbrt(numMolecules)
scale = 5
offset = -0.2
fontSize = 14
plt.rcParams.update({'font.size': fontSize})  
plt.figure(figsize=(1.05*scale, 1*scale))
plt.axhline(y=71.70, linestyle='--', color='grey', label='Experiment')
plt.errorbar(cube_roots, meanValues, yerr=stdErrors, fmt='o', color='black', capsize=5, ecolor='black', label=label)
# Add text labels for the number of samples
for i, txt in enumerate(n_samples):
    plt.text(cube_roots[i]+1.5*offset, (meanValues[i] + stdErrors[i])-5*offset, '$n$='+str(int(txt)), fontsize=int(0.8*fontSize))
# Add text to specify the type of error bars
plt.text(cube_roots[0]+1.5*offset, (np.min(meanValues) - np.max(stdErrors))+40*offset, 'error bar: standard error; $n$: sample size', fontsize=int(0.8*fontSize))
plt.xticks(cube_roots, labels=numMolecules)  # Label ticks with the actual cubic numbers
plt.xlabel('Number (N) of water molecules')
plt.ylabel('Surface tension $\gamma$ (mN/m)')
plt.tick_params(direction="in", axis='both', which='both', top=True, right=True)
plt.xlim(4.5, 10.5)
plt.ylim(40, 100)
plt.legend(frameon=False, ncol=2, loc='upper left')
plt.tight_layout()
plt.savefig('ST_{}.png'.format(simType), dpi=300)
plt.savefig('ST_{}.pdf'.format(simType))
#plt.show()


# ----------------------------
# Plot the surface tension values
# with sliding and accumulated averages
# ----------------------------
plt.rcParams.update({'font.size': fontSize})
fig, axes = plt.subplots(3, 2, figsize=(3 * scale, 1.5 * scale), sharex=True)
# Flatten the 2D axes array for easier iteration
axes = axes.flatten()

for idx, (num, ax) in enumerate(zip(numMolecules, axes)):
    # Plot sliding values
    sliding_filename = "output/ST_Sliding_{}_{}.txt".format(num, simType)
    sliding_data = np.loadtxt(sliding_filename)
    time_sliding, st_values_sliding = sliding_data[:, 0], sliding_data[:, 1]
    ax.plot(time_sliding, st_values_sliding, linewidth=1, label='Sliding Average', color='black', marker='o', markersize=2.5)

    # Plot accumulated values
    accumulated_filename = "output/ST_accumulated_{}_{}.txt".format(num, simType)
    accumulated_data = np.loadtxt(accumulated_filename)
    time_accumulated, accumulated_st = accumulated_data[:, 0], accumulated_data[:, 1]
    ax.plot(time_accumulated, accumulated_st, label='Accumulated Average', color='#57b796')

    # Horizontal line for experiment
    ax.axhline(y=71.70, linestyle='--', color='grey')
    # Add label in the top right corner
    ax.text(0.95, 0.95, 'N={}'.format(num), transform=ax.transAxes, 
            fontsize=fontSize, verticalalignment='top', horizontalalignment='right')

    # Formatting each subplot
    ax.set_ylim(40, 100)
    ax.set_xlim(0, 4000) if simType == 'TIP4P2005' else ax.set_xlim(0, 800)
    ax.tick_params(direction="in", axis='both', which='both', top=True, right=True)
    ax.legend(frameon=False, ncol=1, loc='upper left')
    if idx in [4, 5]:
        ax.set_xlabel('Time (ps)')

# Shared labels for x and y axis
fig.text(0.01, 0.5, 'Surface tension $\gamma$ (mN/m)', ha='center', va='center', rotation='vertical')

plt.tight_layout()
plt.subplots_adjust(hspace=0.1, wspace=0.1)  # Adjust space between plots
plt.savefig('ST_Sliding_Accumulated_{}.png'.format(simType), dpi=300)
plt.savefig('ST_Sliding_Accumulated_{}.pdf'.format(simType))

if simType == 'MB-pol':
    # ----------------------------
    # Plot the surface tension values
    # ----------------------------
    plt.rcParams.update({'font.size': fontSize})
    fig, axes = plt.subplots(3, 1, figsize=(1.7 * scale, 1.5 * scale), sharex=True)
    # Flatten the 2D axes array for easier iteration
    axes = axes.flatten()
    for idx, (num, ax) in enumerate(zip((125, 216, 512), axes)): # Only plot for 125, 216, and 512 molecules
        # Plot sliding values
        sliding_filename = "output/ST_Sliding_{}_{}.txt".format(num, simType)
        sliding_data = np.loadtxt(sliding_filename)
        time_sliding, st_values_sliding = sliding_data[:, 0], sliding_data[:, 1]
        ax.plot(time_sliding, st_values_sliding, linewidth=1, label='Sliding Average', color='black', marker='o', markersize=2.5)
    
        # Plot accumulated values
        accumulated_filename = "output/ST_accumulated_{}_{}.txt".format(num, simType)
        accumulated_data = np.loadtxt(accumulated_filename)
        time_accumulated, accumulated_st = accumulated_data[:, 0], accumulated_data[:, 1]
        ax.plot(time_accumulated, accumulated_st, label='Accumulated Average', color='#57b796')
    
        # Horizontal line for experiment
        ax.axhline(y=71.70, linestyle='--', color='grey')
        # Add label in the top right corner
        ax.text(0.95, 0.95, 'N={}'.format(num), transform=ax.transAxes,
                fontsize=fontSize, verticalalignment='top', horizontalalignment='right')
    
        # Formatting each subplot
        ax.set_ylim(40, 100)
        ax.set_xlim(0, 4000) if simType == 'TIP4P2005' else ax.set_xlim(0, 800)
        ax.tick_params(direction="in", axis='both', which='both', top=True, right=True)
        ax.legend(frameon=False, ncol=1, loc='upper left')
        if idx in [2]:
            ax.set_xlabel('Time (ps)')
    
    # Shared labels for x and y axis
    fig.text(0.015, 0.5, 'Surface tension $\gamma$ (mN/m)', ha='center', va='center', rotation='vertical')
    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)  # Adjust space between plots
    plt.savefig('ST_Sliding_Accumulated_{}_.png'.format(simType), dpi=300)
    plt.savefig('ST_Sliding_Accumulated_{}_.pdf'.format(simType))

# ----------------------------
# Plot the surface tension values 
# with 95% confidence interval
# ----------------------------
meanValues = np.zeros(numMolecules.shape)
stdValues = np.zeros(numMolecules.shape)
stdErrors = np.zeros(numMolecules.shape)
for idx, num in enumerate(numMolecules):
    sliding_filename = "output/ST_Sliding_{}_{}.txt".format(num, simType)
    sliding_data = np.loadtxt(sliding_filename)
    time_sliding, st_values_sliding = sliding_data[:, 0], sliding_data[:, 1]
    meanValues[idx] = np.mean(st_values_sliding)
    stdValues[idx] = np.std(st_values_sliding)
    stdErrors[idx] = stdValues[idx] / np.sqrt(len(st_values_sliding))

confidence_interval = 2 * stdErrors
cube_roots = np.cbrt(numMolecules)
scale = 5
offset = -0.2
fontSize = 14
plt.rcParams.update({'font.size': fontSize})  
plt.figure(figsize=(1.05*scale, 1*scale))
plt.axhline(y=71.70, linestyle='--', color='grey', label='Experiment')
plt.errorbar(cube_roots, meanValues, yerr=confidence_interval, fmt='o', color='black', capsize=5, ecolor='black', label=label)
# Add text labels for the number of samples
for i, txt in enumerate(n_samples):
    plt.text(cube_roots[i]+1.5*offset, (meanValues[i] + confidence_interval[i])-3*offset, '$n$='+str(int(txt)), fontsize=int(0.8*fontSize))
# Add text to specify the type of error bars
plt.text(cube_roots[0]+1.5*offset, (np.min(meanValues) - np.max(confidence_interval))+40*offset, 'error bar: 95% confidence interval; $n$: sample size', fontsize=int(0.8*fontSize))
plt.xticks(cube_roots, labels=numMolecules)  # Label ticks with the actual cubic numbers
plt.xlabel('Number (N) of water molecules')
plt.ylabel('Surface tension $\gamma$ (mN/m)')
plt.tick_params(direction="in", axis='both', which='both', top=True, right=True)
plt.xlim(4.5, 10.5)
plt.ylim(40, 100)
plt.legend(frameon=False, ncol=2, loc='upper left')
plt.tight_layout()
plt.savefig('ST95_{}.png'.format(simType), dpi=300)
plt.savefig('ST95_{}.pdf'.format(simType))
plt.show()

