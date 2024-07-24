#!/usr/bin/env python
import sys 
import json
from ase.io import read
from tqdm import tqdm

def writeOneFrameTo(f, numAtoms, i, time, elements, positions):
    """
    Write one frame to the file f
    """
    assert len(elements) == numAtoms
    assert len(positions) == numAtoms
    f.write("     {}\n".format(numAtoms))
    f.write(" i = {:8d}, time = {:12.3f}\n".format(i, time))
    for j in range(numAtoms):
        f.write("  {}  {:19.10f} {:19.10f} {:19.10f}\n".format(elements[j], positions[j][0], positions[j][1], positions[j][2]))

def extractInfo(logPath):
    with open(logPath, 'r') as file:
        lines = file.readlines()
    sizeX, sizeY, sizeZ = None, None, None
    numAtoms = None

    for idx, line in enumerate(lines):
        if "orthogonal box" in line:
            parts = line.split()
            sizeX = float(parts[-3].replace('(', ''))  # Clean and convert to float
            sizeY = float(parts[-2].replace(')', ''))  # Clean and convert to float
            sizeZ = float(parts[-1].replace(')', ''))  # Clean and convert to float
        if "atoms in group O" in line:
            numAtoms = 3 * int(line.split(' ')[0])  # Extracting the number of atoms
        if "reading angles ..." in line:
            angle_num = lines[idx+1].split(' ')[-2]
            O_num = int(angle_num) # Update the O_num
            print('H2O number: {}'.format(O_num))
            numAtoms = 3 * O_num

    # Ensure that all values were found
    if sizeX is None or sizeY is None or sizeZ is None:
        raise ValueError("The box dimensions could not be found in the log file.")
    if numAtoms is None:
        raise ValueError("The number of atoms could not be found in the log file.")
    return sizeX, sizeY, sizeZ, numAtoms

if __name__ == '__main__':
    if len(sys.argv) == 3: 
        xyzfile = sys.argv[1]
        # Log file log.lammps is in the same directory as the dump.xyz file
        logPath = xyzfile.replace('dump.xyz', 'log.lammps')
        system = sys.argv[2].replace('.xyz', '')
        name = sys.argv[2]
        simType = 'MB-pol'
        # Print xyz file and log file paths
        print("XYZ file: ", xyzfile)
        print("Log file: ", logPath)
    else:
        print('Usage: {} <input.xyz> <output.xyz>'.format(sys.argv[0]))
        exit()
    
    time_step = 0.5 # fs
    dump_every = 1 # dump every dump_freq steps
    ns = 80 # sampling every ns steps
    #oriSimTime = 120 # ps 
    #oriNumFrames = int(oriSimTime * 1000 / time_step / dump_freq) # 24,0000
    
    # Read the the last 100 ps of the trajectory
    # cuttedSimTime = 100 # ps
    cuttedSimTime = 40.0 # 40 ps for testing
    lastNumFrames = int( cuttedSimTime * 1000 / time_step / dump_every) + 1 # 20,0000
   
    # --- For the output JSON file ---
    dt = time_step * dump_every * ns / 1000 # ps
    simTime = cuttedSimTime
    numFrames = int(simTime / dt) + 1 # numFrames after resampling
    sizeX, sizeY, sizeZ, numAtoms = extractInfo(logPath)
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
        "simTime": simTime,
        "numFrames": numFrames
        }

    with open("{}.json".format(system), 'w') as file:
        json.dump(data, file, indent=4)
    print("JSON file {}.json generated successfully.".format(system))
    print("--------------------")

    
    # --- Cut and resample the trajectory ---
    print('Loading xyz file...')
    traj = read(xyzfile, index='-{}:'.format(lastNumFrames))
    
    if len(traj) != lastNumFrames:
        raise ValueError('The number of frames is not equal to the expected number of frames')
    
    numAtoms = len(traj[0])
    elements = traj[0].get_chemical_symbols()
    f = open(sys.argv[2], 'w')
    for i in tqdm(range(0, lastNumFrames, ns)):
        positions = traj[i].get_positions()
        t = i * dt * dump_every
        writeOneFrameTo(f, numAtoms, int((i+1)/ns), t, elements, positions)
    f.close()
    print('File {} has been written'.format(sys.argv[2]))
