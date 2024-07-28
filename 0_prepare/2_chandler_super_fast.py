#!/usr/bin/env python3
import os
from ase.io import read, write
from ase.data import atomic_numbers
from ase.visualize import view
import numpy as np 
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed



def distances2(a, b, box, mic=True):
    """
    Calculate the squared distances between point a and point(s) b in a periodic box.

    Parameters:
    - a: numpy array, shape (3,): Position of the point a.
    - b: numpy array, shape (n, 3): Positions of the points b.
    - box: numpy array, shape (3,): Size of the simulation box.
    - mic: bool: If True, apply the minimum image convention for periodic boundary conditions.

    Returns:
    - squared_distances: numpy array, shape (n,): Squared distances between a and each point in b.
    """
    a = np.array(a)
    b = np.array(b)
    box = np.array(box)
    
    # Compute the distance vector
    delta = b - a

    if mic:
        # Apply minimum image convention
        delta = delta - box * np.round(delta / box)

    # Calculate the squared distance
    squared_distances = np.sum(delta**2, axis=1)
    
    return squared_distances


def ask():
    '''
    Ask the user to input the parameters of the simulation or from a file
    '''
    print('---------------------')
    print('Questions:')
    print('---------------------')
    name=input("What is the name of the position file? (input)\n")
    line=input("What are the crystal parameters? (a b c in Ang)\n")
    line_split=line.split(' ')
    a=float(line_split[0])
    b=float(line_split[1])
    c=float(line_split[2])
    box = np.array([a,b,c])

    whish_size=float(input("What is aproximatively the size of a grid division? (in Ang)\n"))
    n_axis=int(input("What is the normal axis? (0->x, 1->y, 2->z)\n"))
    nb_steps=int(input("What is the total number of steps?\n"))
    ns =int(input("What is ns? (integer: Sample once every ns steps)\n"))
    file_i = open(name,"r")
    nb_atoms=int(file_i.readline()) # The first line is the number of atoms
    file_i.close()
    nb_divx=round(a/whish_size)  # number of division along x
    nb_divy=round(b/whish_size)  #                          y
    nb_divz=round(c/whish_size)  #                          z
    divx=a/nb_divx  # length of each grid along x
    divy=b/nb_divy  #                           y
    divz=c/nb_divz  #   
    print('---------------------')
    print('Answers:')
    print('---------------------')
    print('The name of the position file is: ', name)
    print('The crystal parameters are: a =', a, 'b =', b, 'c =', c)
    print('The size of a grid division is:', whish_size)
    print('The normal axis is:', n_axis)
    print('The total number of steps is:', nb_steps)
    print('The number of ns is:', ns)
    print('The number of atoms is:', nb_atoms)
    print('The number of divisions along x is:', nb_divx)
    print('The number of divisions along y is:', nb_divy)
    print('The number of divisions along z is:', nb_divz)
    print('The length of voxel along x is: {:.4f} Angstrams'.format(divx))
    print('The length of voxel along y is: {:.4f} Angstrams'.format(divy))
    print('The length of voxel along z is: {:.4f} Angstrams'.format(divz))
    return name, a, b, c, box,  whish_size, n_axis, nb_steps, ns, nb_atoms, nb_divx, nb_divy, nb_divz, divx, divy, divz

def recenterAndWrap(atoms, basedOn='O', n_axis=2):
    '''
    Recenter the atoms and wrap them. 
    The atoms are recentered based on the center of mass of the selected atoms.
    The atoms are wrapped based on the periodic boundary conditions.
    Inputs:
    - atoms: Atoms object: The atoms to recenter and wrap.
    - basedOn: str: The element to recenter the atoms based on.
    - n_axis: int: The normal axis.
    Returns:
    - atoms: Atoms object: The recentered and wrapped atoms.
    '''
    # Print cell lengths 
    cell = atoms.get_cell()
    box_lengths = np.linalg.norm(cell, axis=1)
    centerOfBox = box_lengths / 2
    #print('The center of box:', centerOfBox)
    selectedAtoms = atoms[atoms.numbers == atomic_numbers[basedOn]]
    centerOfMass = selectedAtoms.get_center_of_mass()
    #print('The center of mass is:', centerOfMass)

    # Get the center of simulation box
    offset = centerOfBox - centerOfMass

    # Recenter along the normal axis
    atoms.positions[:,n_axis] += offset[n_axis]

    # Warp the positions of the atoms based on the  periodic boundary conditions
    atoms.wrap()

    return atoms

def rhoOAtOrigin(nb_divx, nb_divy, nb_divz, xi=2.4):
    '''
    Calculate the density around a O atom at the origin.
    Inputs:
    - nb_divx: int: The number of divisions along x.
    - nb_divy: int: The number of divisions along y.
    - nb_divz: int: The number of divisions along z.
    - xi: float: The width of the gaussian function.
    Returns:
    - rho: numpy array, shape (nb_divz, nb_divy, nb_divx): The density of the grid points around the O
    '''
    rho = np.zeros((nb_divz,nb_divy,nb_divx))
    #print('The shape of rho is:', rho.shape)
    #print('The shape of rho[0]:', rho[0].shape)

    # Reshape the rho array to a 1D array
    rho = rho.reshape(-1) # Need to be checked

    # Calculate the density around a O atom at a given position
    # Using gaussian function to fill the values of rho around the O atom
    origin = np.array([0,0,0])
    O_pos = origin
    #print('The position of the O atom is:', O_pos)

    # Get all the positions of grid points
    #grids_idx = np.array([[k, j, i] for k in range(nb_divz) for j in range(nb_divy) for i in range(nb_divx)])
    grids_pos = np.array([[i*divx + divx/2, j*divy + divy/2, k*divz + divz/2] for k in range(nb_divz) for j in range(nb_divy) for i in range(nb_divx)])
    #print('The shape of grid_idx is:', grids_idx.shape)
    #print('The shape of grid_points is:', grids_pos.shape)

    # Calculate the squared distances between the O atom and all the grid points
    distances = distances2(O_pos, grids_pos, box)
    #print('The shape of distances is:', distances.shape)


    #xi = 2.4 # Angstrom
    mask = distances < 9*xi**2
    #print('The number of grid points within 3*xi is:', mask.sum())
    #print('The percentage of grid points within 3*xi is: {:.2f}%'.format(mask.sum()/len(mask)*100))


    # Calculate the density of the grid points within 3*xi
    rho[mask] = np.exp(-distances[mask]/(2*xi**2))

    # Reshape back the rho array to 3D
    rho = rho.reshape((nb_divz,nb_divy,nb_divx))

    dim = 3
    normFactor = (2*np.pi*xi**2)**(-dim/2)
    rho = rho * normFactor
    return rho

def calRhoAt(position, rho0, divx, divy, divz):
    '''
    Calculate the density at a given position.
    Inputs:
    - position: numpy array, shape (3,): The position of the atom.
    - rho0: numpy array, shape (nb_divz, nb_divy, nb_divx): The density of the grid points around the O.
    - divx: float: The length of the division along x.
    - divy: float: The length of the division along y.
    - divz: float: The length of the division along z.
    Returns:
    - rho: numpy array, shape (nb_divz, nb_divy, nb_divx):
    '''
    # Move the rho0 array to new position
    origin = np.array([0,0,0])
    offset = position - origin
    shift = np.round(offset/np.array([divx, divy, divz])).astype(int) 
    #print('The error bring by rounding is:', offset/np.array([divx, divy, divz]) - shift)
    #print('The shift is:', shift)

    # Relocate rho0 to the position
    rho = np.roll(rho0, shift, axis=(2,1,0))
    return rho

def visulizeRhoAlongZ(rho, nb_divz, divz):
    '''
    Show the density along the z axis.
    Inputs:
    - rho: numpy array, shape (nb_divz, nb_divy, nb_divx): The density of the grid points around the O.
    - nb_divz: int: The number of divisions along z.
    - divz: float: The length of the division along z.
    '''
    z_length = nb_divz * divz
    slices = np.linspace(0, z_length, 20, endpoint=False)
    # Plot the density, set the same color scale for all the plots
    vmin = rho.min()
    vmax = rho.max()
    fig = plt.figure(figsize=(22, 5))
    gs = gridspec.GridSpec(2, 10)
    axs = [fig.add_subplot(gs[i]) for i in range(20)]
    cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # Adjust the position for the colorbar

    for i, slice in enumerate(slices):
        ax = axs[i]
        im = ax.imshow(rho[int(slice/divz), :, :], vmin=vmin, vmax=vmax, origin='lower')
        ax.axis('off')
        ax.set_title('z = {:.2f} Å'.format(slice))

    # Add the colorbar with the same color scale for all plots
    fig.colorbar(im, cax=cax)
    plt.subplots_adjust(left=0.05, right=0.9, top=0.9, bottom=0.1, wspace=0.1, hspace=0.3)    
    plt.show()

def writeCubeFile(atoms, rho, filename='output_density.cube'):
    '''
    Write the density to a cube file.
    Inputs:
    - atoms: Atoms object: The atoms object.
    - rho: numpy array, shape (nb_divz, nb_divy, nb_divx): The density of atoms.
    - filename: str: The name of the output file.
    '''
    data = np.transpose(rho, (2, 1, 0))
    write(filename, atoms, data=data, format='cube')

def compute_surface(ix, iy, rho, z_positions, density_threshold, nb_divz):
    '''
    Compute the surface of the density profile at a given position.
    Inputs:
    - ix: int: The index along x.
    - iy: int: The index along y.
    - rho: numpy array, shape (nb_divz, nb_divy, nb_divx): The density of the grid points around the O.
    - z_positions: numpy array, shape (nb_divz,): The positions along z.
    - density_threshold: float: The density threshold.
    - nb_divz: int: The number of divisions along z.
    '''
    density_profile = rho[:, iy, ix]
    above_threshold = density_profile >= density_threshold
    surf1_val = surf2_val = np.nan

    if np.any(above_threshold):
        crossing_indices = np.where(above_threshold)[0]
        lower_index = crossing_indices[0]
        upper_index = crossing_indices[-1]

        if lower_index > 0:
            f_lower = interp1d(density_profile[lower_index-1:lower_index+1], 
                               z_positions[lower_index-1:lower_index+1], 
                               kind='linear', bounds_error=False, fill_value=(z_positions[0], z_positions[-1]))
            surf1_val = f_lower(density_threshold)
        else:
            surf1_val = z_positions[lower_index]

        if upper_index < nb_divz - 1:
            f_upper = interp1d(density_profile[upper_index:upper_index+2], 
                               z_positions[upper_index:upper_index+2], 
                               kind='linear', bounds_error=False, fill_value=(z_positions[0], z_positions[-1]))
            surf2_val = f_upper(density_threshold)
        else:
            surf2_val = z_positions[upper_index]
    
    return ix, iy, surf1_val, surf2_val


def obtainSurfaces(rho, nb_divx, nb_divy, nb_divz, divz, density_threshold):
    '''
    Calculate the surfaces of the density profile.
    Inputs:
    - rho: numpy array, shape (nb_divz, nb_divy, nb_divx): The density of the grid points around the O.
    - nb_divx: int: The number of divisions along x.
    - nb_divy: int: The number of divisions along y.
    - nb_divz: int: The number of divisions along z.
    - divz: float: The length of the division along z.
    - density_threshold: float: The density threshold.
    - n_axis: int: The normal axis.
    '''
    z_positions = np.arange(nb_divz) * divz
    surf1 = np.full((nb_divy, nb_divx), np.nan)
    surf2 = np.full((nb_divy, nb_divx), np.nan)

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(compute_surface, ix, iy, rho, z_positions, density_threshold, nb_divz) 
                   for ix in range(nb_divx) for iy in range(nb_divy)]

        for future in futures:
            ix, iy, surf1_val, surf2_val = future.result()
            surf1[iy, ix] = surf1_val
            surf2[iy, ix] = surf2_val

    return surf1, surf2

def visualize_surfaces(surf1, surf2, divx, divy):
    nb_divy, nb_divx = surf1.shape
    x = np.arange(nb_divx) * divx
    y = np.arange(nb_divy) * divy
    X, Y = np.meshgrid(x, y)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the first surface
    ax.plot_surface(X, Y, surf1, color='blue', alpha=0.6, label='Surface 1')
    ax.plot_wireframe(X, Y, surf1, color='blue', alpha=0.3)

    # Plot the second surface
    ax.plot_surface(X, Y, surf2, color='red', alpha=0.6, label='Surface 2')
    ax.plot_wireframe(X, Y, surf2, color='red', alpha=0.3)

    ax.set_xlabel('X-axis (Å)')
    ax.set_ylabel('Y-axis (Å)')
    ax.set_zlabel('Z-axis (Å)')
    ax.set_title('3D Visualization of Density Surfaces')
    plt.legend()
    plt.show()

def writeSurfaceToFile(idx, surf1, surf2, filename='surf.dat'):
    '''
    Write the surfaces to a file.
    Inputs:
    - idx: int, the index of the frame.
    - surf1: numpy array, shape (nb_divy, nb_divx): The first surface.
    - surf2: numpy array, shape (nb_divy, nb_divx): The second surface.
    - filename: str: The name of the output file.
    '''
    nb_divy, nb_divx = surf1.shape
    with open(filename, 'a') as file:
        file.write('i = {}\n'.format(idx))
        for iy in range(nb_divy):
            for ix in range(nb_divx):
                file.write(" {0:5d} {1:5d}{2:12.6f} {3:12.6f}\n".format(iy, ix, surf1[iy, ix], surf2[iy, ix]))

def process_frame(idx, atoms, box, n_axis, rho0, divx, divy, divz, density_threshold):
    atoms.set_cell(box)
    atoms.set_pbc([True, True, True])
    atoms = recenterAndWrap(atoms, basedOn='O', n_axis=n_axis)

    O_atoms = atoms[atoms.numbers == atomic_numbers['O']]
    total_rho = np.zeros_like(rho0)
    for o_pos in O_atoms.positions:
        rho = calRhoAt(o_pos, rho0, divx, divy, divz)
        total_rho += rho
    
    surf1, surf2 = obtainSurfaces(total_rho, nb_divx, nb_divy, nb_divz, divz, density_threshold)

    return idx, total_rho, surf1, surf2

if __name__ == '__main__':
    # Creat a folder called output if it does not exist
    if not os.path.exists('output'): os.makedirs('output')
    name, a, b, c, box, whish_size, n_axis, nb_steps, ns, nb_atoms, nb_divx, nb_divy, nb_divz, divx, divy, divz = ask()
    surfacesFile = 'output/surf_{}.dat'.format(name.split('.')[0])
    os.remove(surfacesFile) if os.path.exists(surfacesFile) else None
    rho0 = rhoOAtOrigin(nb_divx, nb_divy, nb_divz)
    #visulizeRhoAlongZ(rho0, nb_divz, divz)

    # Read the atoms trajectory
    traj = read(name, index=':')
    density_threshold = 0.016  # Example threshold value
    results = []

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_frame, idx, atoms, box, n_axis, rho0, divx, divy, divz, density_threshold)
                   for idx, atoms in enumerate(traj)]

        for future in tqdm(as_completed(futures), total=len(futures)):
            idx, total_rho, surf1, surf2 = future.result()
            results.append((idx, total_rho, surf1, surf2))
            if idx % 10 == 0:
                cubeFile = 'output/{}_{}.cube'.format(name.split('.')[0], idx)
                writeCubeFile(traj[idx], total_rho, filename=cubeFile)

    # Sort results by index to ensure correct order
    results.sort(key=lambda x: x[0])
    for idx, total_rho, surf1, surf2 in results:
        writeSurfaceToFile(idx, surf1, surf2, filename=surfacesFile)

    # ------------------------------------------------------
    # If you want to run the code without multiprocessing,
    # you can use the following codes, which tells the main
    # idea of the code.
    # ------------------------------------------------------
    '''
    for idx, atoms in enumerate(tqdm(traj)):
       #view(atoms)
       atoms.set_cell(box)
       atoms.set_pbc([True,True,True])
       atoms = recenterAndWrap(atoms, basedOn='O', n_axis=n_axis)

       # a) Calculate density
       O_atoms = atoms[atoms.numbers == atomic_numbers['O']]
       #view(O_atoms)
       total_rho = np.zeros_like(rho0)
       for o_pos in O_atoms.positions:
          rho = calRhoAt(o_pos, rho0, divx, divy, divz)
          total_rho = total_rho + rho
       #visulizeRhoAlongZ(total_rho, nb_divz, divz)

       if idx % 10 == 0:
          cubeFile = '{}_{}.cube'.format(name.split('.')[0], idx)
          writeCubeFile(atoms, total_rho, filename=cubeFile)

       # b) Calculate the surfaces
       surf1, surf2 = obtainSurfaces(total_rho, nb_divx, nb_divy, nb_divz, divz)
       #print('The shape of surf1 is:', surf1.shape)
       #print('The shape of surf2 is:', surf2.shape)
       #visualize_surfaces(surf1, surf2, divx, divy)

       writeSurfaceToFile(idx, surf1, surf2, filename=surfacesFile)
       '''
