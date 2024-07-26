#!/usr/bin/env python3
from ase.io import read, write
from ase.visualize import view
import periodictable 
import numpy as np 
import math

#==========
# Constants
#==========
xi = 2.4 # Angstrom

#=============================
#Basic questions (6 Questions)
#=============================
name=input("What is the name of the position file? (input)\n")
line=input("What are the crystal parameters? (a b c in Ang)\n")
line_split=line.split(' ')
a=float(line_split[0])
b=float(line_split[1])
c=float(line_split[2])

whish_size=float(input("What is aproximatively the size of a grid division? (in Ang)\n"))
n_axis=int(input("What is the normal axis? (0->x, 1->y, 2->z)\n"))
nb_steps=int(input("What is the total number of steps?\n"))
ns = int(input("What is ns? (integer: Sample once every ns steps)\n"))
file_i = open(name,"r")
nb_atoms=int(file_i.readline()) # The first line is the number of atoms
file_i.close()
print('-------------------')
print('The name of the position file is: ', name)
print('The crystal parameters are: a =', a, 'b =', b, 'c =', c)
print('The size of a grid division is:', whish_size)
print('The normal axis is:', n_axis)
print('The total number of steps is:', nb_steps)
print('The number of ns is:', ns)
print('The number of atoms is:', nb_atoms)
#END OF QUESTIONS

# Some parameters
nb_divx=round(a/whish_size)  # number of division along x
nb_divy=round(b/whish_size)  #                          y
nb_divz=round(c/whish_size)  #                          z
divx=a/nb_divx  # length of each grid along x
divy=b/nb_divy  #                           y
divz=c/nb_divz  #   

print('The number of divisions along x is:', nb_divx)
print('The number of divisions along y is:', nb_divy)
print('The number of divisions along z is:', nb_divz)
print('The length of each grid along x is:', divx)
print('The length of each grid along y is:', divy)
print('The length of each grid along z is:', divz)

#Calculation of the number of boxes to take into account
#Depends on :
#     1)the orientation of the normal axis
#     2)the ratio between the cell parameters and 3xi
if n_axis == 0:
    range_x=0
    range_y=math.ceil(3*xi/b) # ceil
    range_z=math.ceil(3*xi/c)
    surf1_nd = np.zeros((nb_steps,nb_divy,nb_divz))
    surf2_nd = np.zeros((nb_steps,nb_divy,nb_divz))
elif n_axis ==1:
    range_x=math.ceil(3*xi/a) 
    range_y=0
    range_z=math.ceil(3*xi/c)
    surf1_nd = np.zeros((nb_steps,nb_divx,nb_divz))
    surf2_nd = np.zeros((nb_steps,nb_divx,nb_divz))
else :
    range_x=math.ceil(3*xi/a)
    range_y=math.ceil(3*xi/b)
    range_z=0
    surf1_nd = np.zeros((nb_steps,nb_divx,nb_divy))
    surf2_nd = np.zeros((nb_steps,nb_divx,nb_divy))

print('The range along x is:', range_x)
print('The range along y is:', range_y)
print('The range along z is:', range_z)

# Read teh last frame
atoms = read(name)
# add box size to the trajectory
atoms.set_cell([a,b,c])
atoms.set_pbc([True,True,True])

print('The atoms are:', atoms)
# Select O atoms 
O_atoms = atoms[atoms.numbers == 8]

# View the O atoms
print('The O atoms are:', O_atoms)
view(O_atoms)

# Write atoms to xyz file
atoms.write('LastFrame.xyz')

# Write the guassian cube file
atoms.write('LastFrame.cube')

# Read trajectory
#traj = read(name, index=':')
#print('The length of the trajectory is:', len(traj))
#file_i = open(name,"r")
#name_split=name.split('.')
#name2 = 'surf_{}.dat'.format(name_split[0]) 
#file_o2 = open(name2, 'w') #store the surface arrays into file_o2

# Loop over the trajectory





