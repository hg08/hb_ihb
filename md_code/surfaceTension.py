#!/usr/bin/env python
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the uploaded file
file_path = sys.argv[1] if len(sys.argv) == 2 else 'log.lammps'

time_step = 0.1 # fs
dump_every = 5  # dump every 5 steps

# Check the content of the file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Extract relevant data section from the log file
data_start = None
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
# Ensure L_z was found
if L_z is None:
    raise ValueError("The z-axis length could not be found in the log file.")

# Function to clean and validate data properly
def clean_and_validate(line):
    parts = line.split()
    if len(parts) >= 11:
        try:
            # Convert parts to float to ensure they are all numeric
            parts = [float(part) for part in parts[:11]]
            return parts
        except ValueError:
            return None
    return None

# Use the function to clean and filter the data
cleaned_data_rows = [clean_and_validate(line) for line in lines[data_start + 1:] if clean_and_validate(line) is not None]

# Create the DataFrame with correct columns
column_names = ['Step', 'Temp', 'PotEng', 'TotEng', 'Pxx', 'Pxy', 'Pxz', 'Pyy', 'Pyz', 'Pzz', 'Volume']
data_df = pd.DataFrame(cleaned_data_rows, columns=column_names)

# Plot Temp vs Step
plt.figure(figsize=(14, 8))
plt.plot(data_df['Step']*time_step/1000, data_df['Temp'], label='Temperature')
plt.xlabel('Time (ps)')
plt.ylabel('Temperature (K)')
plt.title('Temperature vs Step (H2O number: {})'.format(O_num))
plt.legend()
plt.grid(True)
plt.savefig('Temp.vs.step.{}_H2O.pdf'.format(O_num))
plt.savefig('Temp.vs.step.{}_H2O.png'.format(O_num))
plt.show()

# Plot PotEng vs Step
plt.figure(figsize=(14, 8))
plt.plot(data_df['Step']*time_step/1000, data_df['PotEng'], label='Potential Energy')
plt.xlabel('Time (ps)')
plt.ylabel('Potential Energy (kcal/mol)')
plt.title('Potential Energy vs Time (H2O number: {})'.format(O_num))
plt.legend()
plt.grid(True)
plt.savefig('PotEng.vs.step.{}_H2O.pdf'.format(O_num))
plt.savefig('PotEng.vs.step.{}_H2O.png'.format(O_num))
plt.show()

# Get rid of unstable values
data_df = data_df[data_df['Step'] > 5000]

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

# Apply sliding window calculation
window_time = 10 # ps
window_size = int(window_time * 1000 / dump_every) 
surface_tension_values = []
steps = []
jump = 20
for i in range(0, len(data_df) - window_size + 1, jump):
    window_slice = data_df.iloc[i:i + window_size]
    surface_tension = calculate_surface_tension(window_slice)
    surface_tension_values.append(surface_tension)
    steps.append(window_slice['Step'].iloc[-1])

times = np.array(steps) * time_step / 1e3 # Time in ps
# Plotting the surface tension values as a function of step
plt.figure(figsize=(14, 8))
plt.plot(times, surface_tension_values, label='Surface Tension (mN/m)')
plt.axhline(y=71.73, color='r', linestyle='--', label='Surface Tension Exp. 300 K: 71.73 mN/m')
plt.xlabel('Time (ps)')
plt.ylabel('Surface Tension (mN/m)')
plt.title('Surface Tension vs Time (H2O number: {}, Sliding Window Size {} ps)'.format(O_num, window_time))
plt.legend()
plt.grid(True)
plt.savefig('ST_Sliding_{}.pdf'.format(O_num))
plt.savefig('ST_Sliding_{}.png'.format(O_num))
plt.show()

# Calculate the mean value changes with time
surface_tension_values = []
for i in range(0, len(data_df) - window_size + 1, jump):
    window_slice = data_df.iloc[0:i + window_size] # From start to now
    surface_tension = calculate_surface_tension(window_slice)
    surface_tension_values.append(surface_tension)

# Plotting the surface tension values as a function of step
plt.figure(figsize=(14, 8))
plt.plot(times, surface_tension_values, label='Surface Tension (mN/m)')
plt.axhline(y=71.73, color='r', linestyle='--', label='Surface Tension Exp. 300 K: 71.73 mN/m')
plt.xlabel('Time (ps)')
plt.ylabel('Surface Tension (mN/m)')
plt.title(' Mean Surface Tension vs Time (H2O number: {})'.format(O_num))
plt.legend()
plt.grid(True)
plt.savefig('Mean_ST_Sliding_{}.pdf'.format(O_num))
plt.savefig('Mean_ST_Sliding_{}.png'.format(O_num))
plt.show()

# Print the value
result = 'Surface tension calculated from simulation: {:.2f} mN/m (Exp. ref. 300 K: 71.73 mN/m)'.format(calculate_surface_tension(data_df))
print(result)

# Open a file in write mode
with open("SurfaceTension.txt", "w") as file:
    # Write the string to the file
    file.write(result)

