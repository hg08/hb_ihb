#!/usr/bin/env python
import sys
import numpy as np
from ase.io import read
from ase.data import atomic_numbers
from tqdm import tqdm

def writeOneFrameTo(f, numAtoms, i, time, elements, positions, box):
    """
    Write one frame to the file f
    """
    boxinfo = "Lattice=\"{} 0.0 0.0 0.0 {} 0.0 0.0 0.0 {}\"".format(box[0], box[1], box[2])
    assert len(elements) == numAtoms
    assert len(positions) == numAtoms
    f.write("     {}\n".format(numAtoms))
    #f.write(" i = {:8d}, time = {:12.3f} {}\n".format(i, time, boxinfo))
    f.write(" i = {:8d}, time = {:12.3f}\n".format(i, time))
    for j in range(numAtoms):
        f.write("  {}  {:19.10f} {:19.10f} {:19.10f}\n".format(elements[j], positions[j][0], positions[j][1], positions[j][2]))

def recenterAndWrap(atoms, box, offset, n_axis=2):
    '''
    Recenter the atoms to the center of the box and wrap
    input: atoms, box, offset, n_axis
    output: atoms
    '''
    atoms.set_cell(box)
    atoms.set_pbc([True, True, True])
    atoms.wrap() # Make sure x, y directions are wrapped correctly
    atoms.positions[:,n_axis] += offset[n_axis]
    atoms.wrap()  
    return atoms

if __name__ == '__main__':
    if len(sys.argv) == 3:
        xyzfile = sys.argv[1]
    else:
        print('Usage: {} <input.xyz> <output.xyz>'.format(sys.argv[0]))
        exit()

    dt = 0.5 # fs
    dump_freq = 1 # dump every dump_freq steps
    ns = 80 # sampling every ns steps
    #oriSimTime = 120 # ps
    #oriNumFrames = int(oriSimTime * 1000 / dt / dump_freq) # 24,0000

    # Read the the last 120 ps of the trajectory
    cuttedSimTime = 120 # ps
    lastNumFrames = int( cuttedSimTime * 1000 / dt / dump_freq) + 1 

    print('Loading xyz file...')
    traj = read(xyzfile, index='-{}:'.format(lastNumFrames))

    if len(traj) != lastNumFrames:
        raise ValueError('The number of frames is not equal to the expected number of frames')

    atoms = traj[-1]
    sizeX = 15.6404
    sizeY = 15.6404
    sizeZ = 31.2808
    box = np.array([sizeX, sizeY, sizeZ]) 
    atoms.set_cell(box)
    atoms.set_pbc([True, True, True])
    atoms.wrap()

    centerOfBox = box / 2
    selectedAtoms = atoms[atoms.numbers == atomic_numbers['O']]
    centerOfMass = selectedAtoms.get_center_of_mass()
    offset = centerOfBox - centerOfMass

    print('Center of mass: {}'.format(centerOfMass))
    print('Offset: {}'.format(offset))
    
    numAtoms = len(traj[0])
    elements = traj[0].get_chemical_symbols()
    f = open(sys.argv[2], 'w')
    for i in tqdm(range(0, lastNumFrames, ns)):
        atoms = recenterAndWrap(traj[i], box, offset)
        positions = atoms.get_positions()
        t = i * dt * dump_freq
        writeOneFrameTo(f, numAtoms, int((i+1)/ns), t, elements, positions, box)
    f.close()
    print('File {} has been written'.format(sys.argv[2]))
