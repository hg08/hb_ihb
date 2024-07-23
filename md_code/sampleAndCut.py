#!/usr/bin/env python
import sys 
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


if __name__ == '__main__':
    if len(sys.argv) == 3: 
        xyzfile = sys.argv[1]
    else:
        print('Usage: {} <input.xyz> <output.xyz>'.format(sys.argv[0]))
        exit()
    
    dt = 0.1 # fs
    dump_freq = 5 # dump every dump_freq steps
    ns = 80 # sampling every ns steps
    #oriSimTime = 120 # ps 
    #oriNumFrames = int(oriSimTime * 1000 / dt / dump_freq) # 24,0000
    
    # Read the the last 100 ps of the trajectory
    cuttedSimTime = 100 # ps
    lastNumFrames = int( cuttedSimTime * 1000 / dt / dump_freq) + 1 # 20,0000
    
    print('Loading xyz file...')
    traj = read(xyzfile, index='-{}:'.format(lastNumFrames))
    
    if len(traj) != lastNumFrames:
        raise ValueError('The number of frames is not equal to the expected number of frames')
    
    numAtoms = len(traj[0])
    elements = traj[0].get_chemical_symbols()
    f = open(sys.argv[2], 'w')
    for i in tqdm(range(0, lastNumFrames, ns)):
        positions = traj[i].get_positions()
        t = i * dt * dump_freq
        writeOneFrameTo(f, numAtoms, int((i+1)/ns), t, elements, positions)
    f.close()
    print('File {} has been written'.format(sys.argv[2]))
