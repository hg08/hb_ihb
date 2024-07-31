#Filename: 1_chandler_fast.py 
# Purpose: 
# 1. To find the two  T*N*N array which define the dynamical surfaces Za0, Zb0 (we only consider Za0, and denote it as Z0 for short.)
# 2. We use numpy to accerate the calculation.           

# To run this function: python3 Filename

import periodictable
import math 
import numpy as np

# To compute running time
import time
import timeit
#For parallel
import multiprocessing
from multiprocessing.dummy import Pool 
import platform

#=========
#CONSTANTS
#=========
atob = 1.88971616463 #Factor to convert from Ang to bohr
dim = 3
xi = 2.4 # Angstrom
xi2_2 = 2*xi**2
xi2_9 = 9*xi**2
norm = (2*math.pi*xi**2)**(-dim/2)
rho_ref = 0.016
max_dist = 1.2 #Maximal length of the O-H bond
max_dist2 = max_dist**2
max_dist_no = 1.4  #Maximal lengh of the O-N bond
max_dist2_no = max_dist_no**2  #Maximal lengh of the O-N 

#===========
#SUBROUTINES
#===========
# Calculate the density associated with one O on a grid point (x,y,z) is atomic coordinate.
# PBC are taken into account
# If distance higher than 3xi we do not take it into account
# We need range_x,..., because we need consider the contribution of Atoms in radius of 3xi on density 
# range_x,...: Nr. of periods of a,...
def density(x,y,z,a,b,c,range_x,range_y,range_z): 
    """The O atom and its image's contribution on the density rho(x,y,z) at point (x,y,z). 
       The main thing: consider the image."""
    d=0
    for i in range(-range_x,range_x+1): 
        for j in range(-range_y,range_y+1): 
            for k in range(-range_z,range_z+1):
                r2=(x-a*i)**2+(y-b*j)**2+(z-c*k)**2 
                #The atoms which are too far away (more than 3*xi) are not taken into account
                if r2 <xi2_9:
                    d += math.exp(-r2/xi2_2) #1 op
    return d 

##This subroutine determines the position of rho_ref on the grid and returns it. 
#def pos_surf(start,end,step,rho,itoc):
#    #First grid point with rho>rho_ref
#    for i in range(start,end,step):
#        if rho[i] >= rho_ref:
#            break
#    #First grid point after the water slab with rho<rho_ref
#    for j in range(i+step,end,step):
#        if rho[j] <= rho_ref :
#            break
#
#    #Linear interpolation of the surface position 
#    #1)Interpolation between i and i-1
#    z1= (rho_ref-rho[i])/(rho[i]-rho[i-step])+i
#    z1= (z1-start)/step*itoc
#    #2)Interpolation between j-1 and j
#    z2= (rho_ref-rho[j])/(rho[j]-rho[j-step])+j
#    z2 =(z2-start)/step*itoc
#    
#    return z1,z2

# Define a function to run the loops in parallel.
def grid(nb_divx, nb_divy, nb_divz):
    """Purpose: Create vertices for a 3D Grid."""
    return [(i, j, k) for i in range(nb_divx) for j in range(nb_divy) for k in range(nb_divz)]

def rho_on_grid(rho_nd,vertice):
    #Sum over all the O, because every O atom contributes to the density at one grid (i,j,k).
    for at in range(nb_atoms):
        if lsymb[at].number > 1 : # Only large atoms are selected
            rho_nd[vertice[0],vertice[1],vertice[2]] += density(lx[at]-(vertice[0]+0.5)*divx, ly[at]-(vertice[1]+0.5)*divy, lz[at]-(vertice[2]+0.5)*divz, a, b, c, range_x, range_y, range_z)
    rho_nd[vertice[0], vertice[1], vertice[2]] *= norm

#This subroutine determines the position of rho_ref on the grid and returns it. 
def pos_surf_nd(x1,x2,rho_nd,itoc):
    #First grid point with rho>rho_ref
    if n_axis == 0:
        for i in range(nb_divx):
            if rho_nd[i,x1,x2] >= rho_ref:
                break
        #First grid point after the water slab with rho<rho_ref
        for ii in range(i+1,nb_divx,1):
            if rho_nd[ii,x1,x2] <= rho_ref :
                break

        #Linear interpolation of the surface position 
        #1)Interpolation between i and i-1
        z1= (rho_ref-rho_nd[i,x1,x2])/(rho_nd[i,x1,x2]-rho_nd[i-1,x1,x2])+i
        z1= z1 * itoc
        #2)Interpolation between j-1 and j
        z2= (rho_ref-rho_nd[ii,x1,x2])/(rho_nd[ii,x1,x2]-rho_nd[ii-1,x1,x2])+ii
        z2 =z2 * itoc
    elif n_axis == 1:
        for j in range(nb_divy):
            if rho_nd[x1,j,x2] >= rho_ref:
                break
        #First grid point after the water slab with rho<rho_ref
        for jj in range(j+1,nb_divy,1):
            if rho_nd[x1,j,x2] <= rho_ref :
                break

        #Linear interpolation of the surface position 
        #1)Interpolation between j and j-1
        z1= (rho_ref-rho_nd[x1,j,x2])/(rho_nd[x1,j,x2]-rho_nd[x1,j-1,x2])+j
        z1= z1 * itoc
        #2)Interpolation between jj-1 and jj
        z2= (rho_ref-rho_nd[x1,jj,x2])/(rho_nd[x1,jj,x2]-rho_nd[x1,jj-1,x2])+jj
        z2 =z2 * itoc
    else:
        for k in range(nb_divz):
            if rho_nd[x1,x2,k] >= rho_ref:
                break
        #First grid point after the water slab with rho<rho_ref
        for kk in range(k+1,nb_divz,1):
            if rho_nd[x1,x2,kk] <= rho_ref :
                break

        #Linear interpolation of the surface position 
        #1)Interpolation between i and i-1
        z1= (rho_ref-rho_nd[x1,x2,k])/(rho_nd[x1,x2,k]-rho_nd[x1,x2,k-1])+k
        z1= z1 * itoc
        #2)Interpolation between j-1 and j
        z2= (rho_ref-rho_nd[x1,x2,kk])/(rho_nd[x1,x2,kk]-rho_nd[x1,x2,kk-1])+kk
        z2 =z2 * itoc
    
    return z1,z2

def identify_os():
    # Get the name of the operating system
    os_name = platform.system()
    if os_name == "Darwin":
        # This block will execute for macOS
        # print("This is a macOS system.")
        # Example: Set multiprocessing start method to 'fork'
        import multiprocessing as mp
        mp.set_start_method('fork', force=True)
    elif os_name == "Linux":
        # This block will execute for Linux
        # print("This is a Linux system.")
        # Example: Continue using the default start method ('fork' on Linux)
        pass
    else:
        print("This is not a macOS or Linux system.")

#=============================
#Basic questions (6 Questions)
#=============================
name=input("What is the name of the position file? (input)\n")
file_i = open(name,"r")
line=input("What are the crystal parameters? (a b c in Ang)\n")
line_split=line.split(' ')
a=float(line_split[0])
b=float(line_split[1])
c=float(line_split[2])

whish_size=float(input("What is aproximatively the size of a grid division? (in Ang)\n"))
n_axis=int(input("What is the normal axis? (0->x, 1->y, 2->z)\n"))
nb_steps=int(input("What is the total number of steps?\n"))
ns = int(input("What is ns? (integer: Sample once every ns steps)\n"))

name_split=name.split('.')

name2 = 'surf_{}.dat'.format(name_split[0]) 

file_o2 = open(name2, 'w') #store the surface arrays into file_o2
#END OF QUESTIONS

#================================================
#Reading the main information
nb_atoms=int(file_i.readline()) #Number of atoms (the first line is the number of atoms)

# Some parameters
nb_divx=round(a/whish_size)  # number of division along x
nb_divy=round(b/whish_size)  #                          y
nb_divz=round(c/whish_size)  #                          z
divx=a/nb_divx  # length of each grid along x
divy=b/nb_divy  #                           y
divz=c/nb_divz  #                           z

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
    
#========================================
#nb_atoms lines about the atomic position
#atomic_number   random_float    x   y  z
#========================================
#Skip the commentary line of the xyz file
file_i.readline()

#1)Read all the atomic positions (and symbol)
#and calculate the position of the mass center
lsymb = [None]*nb_atoms # use a list of length nb_atoms to store symb of atoms
lx = [None]*nb_atoms
ly = [None]*nb_atoms
lz = [None]*nb_atoms
lsymb_nd = np.empty((nb_steps,nb_atoms))
lx_nd = np.zeros((nb_steps,nb_atoms))
ly_nd = np.zeros((nb_steps,nb_atoms))
lz_nd = np.zeros((nb_steps,nb_atoms))
sum_mass_nd=np.zeros(nb_steps)
center_pos_nd=np.zeros(nb_steps)
sum_mass=0
center_pos=0

# Go to start of file
file_i.seek(0)

# n_temp is the accelerate rate, with help of n_temp, when we calculated a surface pairs, then we can skip n_temp steps. 
# ns is the same parameter as ns used in the code 'ihb_sulpizi'
#n_temp = 10 * ns  
n_temp = 1 * ns  

# Calculate the surfaces every n_temp time steps
for s in range(0,nb_steps,ns):
    if (s%n_temp == 0) :
        #skip 2 lines
        file_i.readline()
        file_i.readline()
        
        #file_o2.write("i = {0:5d}\n".format(s))    
        file_o2.write("i = {}\n".format(s))    

        name0="".join([name_split[0],str(s),'.cube'])
        file_o = open(name0, 'w') #the cube file as output

        #================================================
        #The first 2 lines of the cube file are written
        file_o.write("cube file for determining the isosurface of a system.\n")
        file_o.write("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")

        #================================================
        #The 3rd line of the cube file are written
        #{0:5d} means: the data number 0 is an integer (d) and will be printed with 5 digits
        #{1:12.6f} means: the data number 1 is float (f) and will be printed with 12 digits whom 6 will be after the .
        file_o.write("{0:5d} {1:12.6f} {2:12.6f} {3:12.6f}\n".format(nb_atoms,0,0,0))

        #=============================
        #the 4th line of the cube file
        #Definition of the grid
        #=============================
        file_o.write("{0:5d} {1:12.6f} {2:12.6f} {3:12.6f}\n".format(nb_divx,divx*atob,0        ,0))
        file_o.write("{0:5d} {1:12.6f} {2:12.6f} {3:12.6f}\n".format(nb_divy,0        ,divy*atob,0))
        file_o.write("{0:5d} {1:12.6f} {2:12.6f} {3:12.6f}\n".format(nb_divz,0        ,0        ,divz*atob))
            
        for i in range(nb_atoms):
            symb,x,y,z=file_i.readline().split()
            #Record the values inside a list with good format
            lsymb[i]=periodictable.elements.symbol(symb) #lsymb is an object
            lx[i]=float(x)
            ly[i]=float(y)
            lz[i]=float(z)
        
            #Calculate the mass center along the normal axis
            if n_axis==0:
                center_pos += lx[i]*lsymb[i].mass
            elif n_axis==1:
                center_pos += ly[i]*lsymb[i].mass
            else: 
                center_pos += lz[i]*lsymb[i].mass
            sum_mass +=lsymb[i].mass 

        #For each step s, we calculate the center_pos
        center_pos/=sum_mass

        #2)Recenter (along the normal axis only) the atoms according to the mass center
        #if n_axis==0:
        #    for i in range(nb_atoms):
        #        lx[i]=lx[i]-center_pos + a/2
        #elif n_axis ==1:
        #    for i in range(nb_atoms):
        #        ly[i]=ly[i]-center_pos + b/2
        #else:
        #    for i in range(nb_atoms):
        #        # lz_new - lz_old = New_z_center_of_mass - z_center_of_mass (New_z_center_of_mass is selected as c/2, the center of the cell along z-axis)
        #        lz[i]=lz[i]-center_pos + c/2 

        #3)Write the data into the cube file
        for i in range(nb_atoms):
            file_o.write(" {0:5d} {1:12.6f} {2:12.6f} {3:12.6f} {4:12.6f}\n".format(lsymb[i].number,0,lx[i]*atob,ly[i]*atob,lz[i]*atob))

        #===============
        #Volumetric data
        #===============

        inc=0
        rho=[0]*nb_divx*nb_divy*nb_divz
        #Loop for the grid: this loop i-j-k has to be paralleled 
        for i in range (nb_divx): #The program starts with this line
            for j in range(nb_divy):
                for k in range(nb_divz):
                    #Sum over all the O, because every O atom contributes to the density at each grid (i,j,k).
                    for at in range(nb_atoms):
                        if lsymb[at].number == 8: # Only the O are selected
                            rho[inc]+=density(lx[at]-(i+0.5)*divx, ly[at]-(j+0.5)*divy, lz[at]-(k+0.5)*divz, a, b, c, range_x, range_y, range_z)
                    rho[inc] *= norm
                    file_o.write(" {0:12.6f}".format(rho[inc]))
                    inc+=1
                    if k % 6 == 5:
                        file_o.write("\n")
                if k%6 != 5:
                    file_o.write("\n")

        #==============
        # parallel code
        #==============
        
        # First, genereate the 3D grid
        vertices = grid(nb_divx,nb_divy,nb_divz)
        # Second, prepare a 3D array: rho_nd
        rho_nd = np.zeros((nb_divx,nb_divy,nb_divz))
        
        if __name__ == '__main__':
            # Run the function to identify the OS
            identify_os()

            p = multiprocessing.Pool(processes = multiprocessing.cpu_count()-1)

            start = time.time()

            #This step solve the problem of 'How to pass elements of vertices to rho_on_grid'.
            for vertice in vertices:
                p.apply_async(rho_on_grid(rho_nd,vertice),[vertice])

            p.close()
            p.join()
            print("Complete")
            end = time.time()
            print('Total time (s)= ' + str(end-start))
        
            # Check the rho_nd
            #print(grid(nb_divx,nb_divy,nb_divz))
            #print("rho from parallel code :{}".format(rho_nd))
            #print("The first element of rho_nd:{}".format(rho_nd[10,10,10]))

        #====================
        #Determine isosurface
        #====================
        if n_axis==0:
            surf1=[0]*nb_divy*nb_divz # Initialization of the list associated with the perpendicular plane; 
                                      # Usually, we can express surf1 at a certain time as an array, but here we choose list to store surf1.
            surf2=[0]*nb_divy*nb_divz
            # For each point of this plane, the point where the reference density is reached is calculated (and stored inside surf1 surf2)
            for j in range(nb_divy):
                for k in range(nb_divz):
                    start_inc= k + nb_divz*j # to generate an index
                    step_inc = nb_divz*nb_divy  
                    end_inc  = start_inc+step_inc*nb_divx
                    surf1_nd[s,j,k], surf2_nd[s,j,k] = pos_surf_nd(j, k, rho_nd, divx) #This function calculates the position of the reference density
                    file_o2.write(" {0:5d} {1:5d}{2:12.6f} {3:12.6f}\n".format(j,k,surf1_nd[s,j,k], surf2_nd[s,j,k])) 
                    #print(start_inc, step_inc, end_inc)
                    #print(surf1[inc])
        elif n_axis==1:
            surf1=[0]*nb_divx*nb_divz
            surf2=[0]*nb_divx*nb_divz
            for i in range(nb_divx):
                for k in range(nb_divz):
                    start_inc= k+nb_divz*(nb_divy*i)
                    step_inc = nb_divz
                    end_inc  = start_inc+step_inc*nb_divy
                    surf1_nd[s,i,k],surf2_nd[s,i,k] = pos_surf_nd(i, k, rho_nd, divy)
                    file_o2.write(" {0:5d} {1:5d} {2:12.6f} {3:12.6f}\n".format(i,k,surf1_nd[s,i,k], surf2_nd[s,i,k])) 
        else:
            surf1=[0]*nb_divx*nb_divy
            surf2=[0]*nb_divx*nb_divy
            for i in range(nb_divx):
                for j in range(nb_divy):
                    start_inc= nb_divz*(j+nb_divy*i)
                    step_inc = 1
                    end_inc  = start_inc+step_inc*nb_divz
                    surf1_nd[s,i,j], surf2_nd[s,i,j] = pos_surf_nd(i, j, rho_nd, divz)
                    file_o2.write("{0:5d} {1:5d} {2:12.6f} {3:12.6f}\n".format(i,j,surf1_nd[s,i,j], surf2_nd[s,i,j])) 
        print("surf1_nd:{}".format(surf1_nd))
        print("surf2_nd:{}".format(surf2_nd))
    else:
        #=======================
        # Print the index number 
        #=======================
        #file_o2.write("i = {0:5d}\n".format(s))    
        file_o2.write("i = {}\n".format(s))    

        #====================
        #Determine isosurface
        #====================
        if n_axis==0:
            surf1=[0]*nb_divy*nb_divz # Initialization of the list associated with the perpendicular plane; 
                                      # Usually, we can express surf1 at a certain time as an array, but here Remi and Gang choose list to store surf1.
            surf2=[0]*nb_divy*nb_divz
            # For each point of this plane, the point where the reference density is reached is calculated (and stored inside surf1 surf2)
            for j in range(nb_divy):
                for k in range(nb_divz):
                    start_inc= k + nb_divz*j # to generate an index
                    step_inc = nb_divz*nb_divy  
                    end_inc  = start_inc+step_inc*nb_divx
                    surf1_nd[s,j,k], surf2_nd[s,j,k] = pos_surf_nd(j, k, rho_nd, divx) #This function calculates the position of the reference density
                    file_o2.write(" {0:5d} {1:5d}{2:12.6f} {3:12.6f}\n".format(j,k,surf1_nd[s,j,k], surf2_nd[s,j,k])) 
        elif n_axis==1:
            surf1=[0]*nb_divx*nb_divz
            surf2=[0]*nb_divx*nb_divz
            for i in range(nb_divx):
                for k in range(nb_divz):
                    start_inc= k+nb_divz*(nb_divy*i)
                    step_inc = nb_divz
                    end_inc  = start_inc+step_inc*nb_divy
                    surf1_nd[s,i,k],surf2_nd[s,i,k] = pos_surf_nd(i, k, rho_nd, divy)
                    file_o2.write(" {0:5d} {1:5d} {2:12.6f} {3:12.6f}\n".format(i,k,surf1_nd[s,i,k], surf2_nd[s,i,k])) 
        else:
            surf1=[0]*nb_divx*nb_divy
            surf2=[0]*nb_divx*nb_divy
            for i in range(nb_divx):
                for j in range(nb_divy):
                    start_inc= nb_divz*(j+nb_divy*i)
                    step_inc = 1
                    end_inc  = start_inc+step_inc*nb_divz
                    surf1_nd[s,i,j], surf2_nd[s,i,j] = pos_surf_nd(i, j, rho_nd, divz)
                    file_o2.write("{0:5d} {1:5d} {2:12.6f} {3:12.6f}\n".format(i,j,surf1_nd[s,i,j], surf2_nd[s,i,j])) 
        #print("surf1_nd:{}".format(surf1_nd))
        #print("surf2_nd:{}".format(surf2_nd))
        
