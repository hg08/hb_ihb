#!/usr/bin/env python3
import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# Read the uploaded file
file_path = sys.argv[1] if len(sys.argv) == 2 else 'log.lammps'

# Create a output folder
if not os.path.exists('output'):
    os.makedirs('output')

time_step = 0.5 # fs
#dump_every = 1  # dump every 5 steps

# Check the content of the file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Extract relevant data section from the log file
data_start = None
simulationType = ''
for idx, line in enumerate(lines):
    if "Step" in line:
        data_start = idx
        break
    if "orthogonal box" in line:
        print('Simulation box:{}'.format(line[:-1]))
        parts = line.split()
        L_z = float(parts[-1][:-1])  # Extract the z size
    if "atoms in group O" in line:
        O_num = line.split(' ')[0]
        print('H2O number: {}'.format(O_num))
    if "reading angles ..." in line:
        angle_num = lines[idx+1].split(' ')[-2]
        O_num = int(angle_num) # Update the O_num
        print('H2O number: {}'.format(O_num))
    if "### TIP4P Potential Parameters ##" in line:
        simulationType = 'TIP4P2005'
        print('Simulation Type: {}'.format(simulationType))
    if "MBX: A many-body" in line:
        simulationType = 'MB-pol'
        print('Simulation Type: {}'.format(simulationType))
    if "thermo_freq" in line and "equal" in line:
        thermo_every = int(line.split(' ')[-1])
        thermo_every = thermo_every
        print('Thermo every: {}'.format(thermo_every))
    if "dump_freq" in line and "equal" in line:
        dump_every = int(line.split(' ')[-1])
        dump_every = dump_every   
        print('Dump every: {}'.format(dump_every))

# Ensure L_z was found
if L_z is None:
    raise ValueError("The z-axis length could not be found in the log file.")

# Function to clean and validate data properly
def clean_and_validate(line):
    parts = line.split()
    if len(parts) >= 16:
        try:
            # Convert parts to float to ensure they are all numeric
            parts = [float(part) for part in parts[:16]]
            return parts
        except ValueError:
            return None
    return None

# Use the function to clean and filter the data
cleaned_data_rows = [clean_and_validate(line) for line in lines[data_start + 1:] if clean_and_validate(line) is not None]

# Create the DataFrame with correct columns
column_names = ['Step', 'Time',  'Temp', 'TotEng', 'KinEng', 'PotEng', 'Enthalpy', 'Density', 'Lx', 'Ly', 'Lz', 'Volume', 'Pxx', 'Pyy', 'Pzz', 'Press'] # For MB-pol log file
data_df = pd.DataFrame(cleaned_data_rows, columns=column_names)

# Plot Temp vs Step
plt.figure(figsize=(14, 8))
plt.plot(data_df['Time']/1000, data_df['Temp'], label='Temperature')
plt.xlabel('Time (ps)')
plt.ylabel('Temperature (K)')
plt.title('Temperature vs Step (H2O number: {})'.format(O_num))
plt.legend()
plt.grid(True)
plt.savefig('output/Temp.vs.step.{}_H2O_mbpol.pdf'.format(O_num))
plt.savefig('output/Temp.vs.step.{}_H2O_mbpol.png'.format(O_num))
#plt.show()

# Plot PotEng vs Step
plt.figure(figsize=(14, 8))
plt.plot(data_df['Time']/1000, data_df['PotEng'], label='Potential Energy')
plt.xlabel('Time (ps)')
plt.ylabel('Potential Energy (kcal/mol)')
plt.title('Potential Energy vs Time (H2O number: {})'.format(O_num))
plt.legend()
plt.grid(True)
plt.savefig('output/PotEng.vs.step.{}_H2O_{}.pdf'.format(O_num, simulationType))
plt.savefig('output/PotEng.vs.step.{}_H2O_{}.png'.format(O_num, simulationType))
#plt.show()

# Get rid of the unstable values at the beginning
data_df = data_df[data_df['Time'] > 20000]

def calculate_surface_tension(data_slice):
    '''
    Calculate surface tension in a time sliding window
    '''
    Pxx_avg = data_slice['Pxx'].mean()
    Pyy_avg = data_slice['Pyy'].mean()
    Pzz_avg = data_slice['Pzz'].mean()
    
    # Convert the pressure components from atm to Pa
    atm_to_Pa = 101325  # 1 atm = 101325 Pa
    Pxx_avg_Pa = Pxx_avg * atm_to_Pa
    Pyy_avg_Pa = Pyy_avg * atm_to_Pa
    Pzz_avg_Pa = Pzz_avg * atm_to_Pa
    
    # Define the length of the simulation cell in the longest direction (along z axis)
    angstrom_to_m = 1e-10  # 1 Å = 1e-10 meters
    L_z_m = L_z * angstrom_to_m
    
    # Calculate the surface tension in N/m
    surface_tension_N_per_m = 0.5 * L_z_m * (Pzz_avg_Pa - 0.5 * (Pxx_avg_Pa + Pyy_avg_Pa))
    
    # Convert the surface tension from N/m to mN/m
    surface_tension_mN_per_m = surface_tension_N_per_m * 1000  # 1 N/m = 1000 mN/m
    
    return surface_tension_mN_per_m

# --------------------------------
# Apply sliding window calculation
# --------------------------------
window_time = 40 # ps
window_size = int(window_time * 1000 / (time_step * thermo_every))  
surface_tension_values = []
steps = []
jump = window_size
for i in range(0, len(data_df) - window_size + 1, jump):
    window_slice = data_df.iloc[i:i + window_size]
    surface_tension = calculate_surface_tension(window_slice)
    surface_tension_values.append(surface_tension)
    steps.append(window_slice['Step'].iloc[-1])

surface_tension_values = np.array(surface_tension_values)
mean_surface_tension = surface_tension_values.mean()
std_surface_tension = surface_tension_values.std()
n = len(surface_tension_values)
std_error = std_surface_tension / np.sqrt(n)
print('Mean Surface Tension: {:.2f} mN/m'.format(mean_surface_tension))
print('Standard Deviation: {:.2f} mN/m'.format(std_surface_tension))
print('Standard Error: {:.2f} mN/m'.format(std_error))
print('Number of data points: {}'.format(n))
# Save mean and std values to a text file using numpy
np.savetxt('output/ST_mean_std_{}_{}.txt'.format(O_num, simulationType), np.array([[mean_surface_tension, std_surface_tension, std_error, n]]), fmt='%.2f')

times = np.array(steps) * time_step / 1e3 # Time in ps
# Plotting the surface tension values as a function of step
plt.figure(figsize=(14, 8))
plt.plot(times, surface_tension_values, label='Surface Tension (mN/m)')
np.savetxt('output/ST_Sliding_{}_{}.txt'.format(O_num, simulationType), np.array([times, surface_tension_values]).T, fmt='%.2f')
plt.axhline(y=71.73, color='r', linestyle='--', label='Surface Tension Exp. 300 K: 71.73 mN/m')
plt.xlabel('Time (ps)')
plt.ylabel('Surface Tension (mN/m)')
plt.title('Surface Tension vs Time (H2O number: {}, Sliding Window Size {} ps)'.format(O_num, window_time))
plt.legend()
plt.grid(True)
plt.savefig('output/ST_Sliding_{}_{}.pdf'.format(O_num, simulationType))
plt.savefig('output/ST_Sliding_{}_{}.png'.format(O_num, simulationType))
#plt.show()

# --------------------------------
# Calculate the acummulated mean 
# value changes with time
# --------------------------------
surface_tension_values = []
jump = 1
steps = []
for i in tqdm(range(0, len(data_df) - window_size + 1, jump)):
    window_slice = data_df.iloc[0:i + window_size] # From start to now
    surface_tension = calculate_surface_tension(window_slice)
    surface_tension_values.append(surface_tension)
    steps.append(window_slice['Step'].iloc[-1])


times = np.array(steps) * time_step / 1e3 # Time in ps
# Plotting the surface tension values as a function of step
plt.figure(figsize=(14, 8))
plt.plot(times, surface_tension_values, label='Surface Tension (mN/m)')
np.savetxt('output/ST_accumulated_{}_{}.txt'.format(O_num, simulationType), np.array([times, surface_tension_values]).T, fmt='%.2f')
plt.axhline(y=71.73, color='r', linestyle='--', label='Surface Tension Exp. 300 K: 71.73 mN/m')
plt.xlabel('Time (ps)')
plt.ylabel('Surface Tension (mN/m)')
plt.title(' Mean Surface Tension vs Time (H2O number: {})'.format(O_num))
plt.legend()
plt.grid(True)
plt.savefig('output/Mean_ST_Sliding_{}_{}.pdf'.format(O_num, simulationType))
plt.savefig('output/Mean_ST_Sliding_{}_{}.png'.format(O_num, simulationType))
#plt.show()

# Print the value
result = 'Surface tension calculated from simulation: {:.2f} mN/m (Exp. ref. 300 K: 71.73 mN/m)'.format(calculate_surface_tension(data_df))
print(result)

# Open a file in write mode
with open("output/SurfaceTension_{}_{}.txt".format(O_num, simulationType), "w") as file:
    # Write the string to the file
    file.write(result)

