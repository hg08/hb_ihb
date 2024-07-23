#!/usr/bin/env python3
import sys, json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def extractInfo(logPath):
    with open(logPath, 'r') as file:
        lines = file.readlines()
    sizeX, sizeY, sizeZ = None, None, None
    numAtoms = None

    for line in lines:
        if "orthogonal box" in line:
            parts = line.split()
            sizeX = float(parts[-3].replace('(', ''))  # Clean and convert to float
            sizeY = float(parts[-2].replace(')', ''))  # Clean and convert to float
            sizeZ = float(parts[-1].replace(')', ''))  # Clean and convert to float
        if "atoms in group O" in line:
            numAtoms = 3 * int(line.split(' ')[0])  # Extracting the number of atoms

    # Ensure that all values were found
    if sizeX is None or sizeY is None or sizeZ is None:
        raise ValueError("The box dimensions could not be found in the log file.")
    if numAtoms is None:
        raise ValueError("The number of atoms could not be found in the log file.")
    # Only keep 4 decimal places, store as string,  add extra 0s if needed
    sizeX = "{:.4f}".format(sizeX)
    sizeY = "{:.4f}".format(sizeY)
    sizeZ = "{:.4f}".format(sizeZ)
    return sizeX, sizeY, sizeZ, numAtoms

if __name__ == "__main__":
    time_step = 0.1 # fs
    dump_every = 5
    ns = 80 # sample every 
    numWaterList = [125, 216, 343, 512, 729, 1000]
    for numWater in numWaterList:
        system = "{}w-TIP4P2005".format(numWater)
        name = "{}.xyz".format(system)
        simType = "MD"
        dt =  time_step * dump_every * ns / 1000 # in ps
        simTime = 100.0 # in ps
        numFrames = int(1 + simTime / dt)
        logPath = "AirWaterInterface{}/run_TIP4P2005_300_hdf5/run_TIP4P2005_300/log.lammps".format(numWater)
        sizeX, sizeY, sizeZ, numAtoms = extractInfo(logPath)

        # Print all the above information
        print("--------------------")
        print("System: {}".format(system))
        print("Name: {}".format(name))
        print("Type: {}".format(simType))
        print("Time step: {} ps".format(dt))
        print("Simulation time: {} ps".format(simTime))
        print("Number of frames: {}".format(numFrames))
        print("Box dimensions: {} x {} x {}".format(sizeX, sizeY, sizeZ))
        print("Number of atoms: {}".format(numAtoms))

        # Generate the JSON file
        data = {
            "name": name,
            "system": system,
            "simType": simType,
            "dt": dt,
            "sizeX": sizeX,
            "sizeY": sizeY,
            "sizeZ": sizeZ,
            "numAtoms": numAtoms,
            "simTime": "{:.2f}".format(simTime),
            "numFrames": numFrames
            }

        with open("{}.json".format(system), 'w') as file:
            json.dump(data, file, indent=4)
        print("JSON file {}.json generated successfully.".format(system))

        traj = "AirWaterInterface{}/run_TIP4P2005_300_hdf5/run_TIP4P2005_300/{}".format(numWater, name)
        
        # Copy to current directory
        os.system("cp {} .".format(traj))




