#/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_data(filename):
    # Manually specify the headers based on your file structure
    headers = ['Step', 'Time', 'Kin', 'Temp', 'Pot', 'Cons', 'UsedTime']

    # Load text file using numpy, skip the first row
    data = np.loadtxt(filename, skiprows=1)

    # Extracting necessary columns
    time = data[:, 1] * 0.001 # Convert fs to ps
    kin_energy = data[:, 2]
    temp = data[:, 3]
    pot_energy = data[:, 4]

    # Plot Time vs. Kinetic Energy, Temperature, Potential Energy in subplots
    plt.figure(figsize=(12, 8))
    plt.subplot(3, 1, 1)
    plt.plot(time, kin_energy, 'r-')
    plt.xlabel('Time (ps)')
    plt.ylabel('Kinetic Energy (Kcal/mol)')
    plt.title('Time vs. Kinetic Energy')

    plt.subplot(3, 1, 2)
    plt.plot(time, temp, 'g-')
    plt.xlabel('Time (ps)')
    plt.ylabel('Temperature (K)')
    plt.title('Time vs. Temperature')
        
    plt.subplot(3, 1, 3)
    plt.plot(time, pot_energy, 'b-')
    plt.xlabel('Time (ps)')
    plt.ylabel('Potential Energy (Kcal/mol)')
    plt.title('Time vs. Potential Energy')
    
    plt.tight_layout()
    plt.show()
    
    plt.savefig('energies.png')



if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
    else:
        filename = sys.argv[1]
        plot_data(filename)

