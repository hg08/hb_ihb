MODULE module_ihb
  IMPLICIT NONE
  ! Purpose: 
  ! I create the module 'module_ihb', putting all the types, functions and subroutines (in interface) into a single module,
  ! instead of write multiple modules. Because if I use those multiple modules, it's hard to compile succesfully the makefile 
  ! since the module A may use module B. I hope this problem can be solved if I merge all the modules into one, and then 
  ! write subroutines in seperate files.

  integer,parameter,PUBLIC :: rk=8                                     ! For setting the type of real variables 
  real(kind=rk),parameter :: rate=0.8d0                                ! For cutting off autocorrelation functions
  integer,parameter,public :: num_coord_max=20                         ! For setting size of the neighbors_indices 
  integer,parameter,public :: num_coord_mol_max=20                     ! For setting size of the wat_neighbors_indices 
  REAL(kind=rk),PARAMETER,PUBLIC :: alpha = 1.d0                       ! For rescaling the radius of spheres (atoms)
  REAL(kind=rk),PARAMETER,PUBLIC :: r_OH_LU2 = 1.44d0                  ! Least Upper (LU) bound of the squared distance between O and H, since r(OH) < 1.2 A
  REAL(kind=rk),PARAMETER,PUBLIC :: r_ON_LU2 = 1.82d0                  ! LU bound of the squared distance between O and N, since rNO(Nitrate) < 1.35 A
  REAL(kind=rk),PARAMETER,PUBLIC :: r_Li_OW_LU2 = 9.00d0               ! Least Upper (LU) bound of the squared distance between Li and OW, since r(Li-OW) < 3.0 A
  REAL(kind=rk),PARAMETER,PUBLIC :: r_Na_OW_LU2 = 10.24d0              ! Least Upper (LU) bound of the squared distance between Na and OW, since r(Na-OW) < 3.2 A
  REAL(kind=rk),PARAMETER,PUBLIC :: r_K_OW_LU2 = 13.96d0               ! Least Upper (LU) bound of the squared distance between K and OW, since r(K-OW) < 3.6 A
  REAL(kind=rk),PARAMETER,PUBLIC :: r_I_OW_LU2 = 20.25d0               ! Least Upper (LU) bound of the squared distance between I and OW, since r(I-OW) < 4.5 A
  REAL(kind=rk),PARAMETER,PUBLIC :: r_NN_OW_LU2 = 25.00d0              ! LU bound of the squared distance between NN and OW, since r_{NN-OW}(Nitrate) < 5 A (TO CHECK)
  REAL(kind=rk),PARAMETER,PUBLIC :: r_ON_OW_LU2 = 17.64d0              ! LU bound of the squared distance between ON and OW, since r_{ON-OW}(Nitrate) < 4.2 A 
  REAL(kind=rk),PARAMETER,PUBLIC :: r_OW_OW_LU2 = 12.25d0              ! LU bound of the squared distance between OW and OW, since r_{OW-OW}(Nitrate) < 3.5 A 
  REAL(KIND=rk), PARAMETER,PUBLIC :: whish_size=0.5d0                  ! uinit:Angstrom; For griding

  !Define the derived data type for atom
  TYPE :: atom
    CHARACTER(LEN=2) :: atom_name
    INTEGER :: atom_id 
    INTEGER :: host_id                                                 ! For O atom in water, host_id = atom_id
    REAL(kind=8) :: mass
    REAL(kind=8), DIMENSION(3) :: coord 
  END TYPE atom

  !Define the derived data type for sphere/atom
  ! Ref: Busa2005's package, ARVO
  TYPE :: sphere
    CHARACTER(LEN=2) :: name
    INTEGER :: id 
    REAL(kind=8) :: mass
    REAL(kind=8) :: radius ! the radius of the sphere 
    REAL(kind=8), DIMENSION(3) :: coord  ! x,y,z: where (x,y,z) is the coordinates of the center of the sphere
  END TYPE sphere

  ! Derived data type to stroe array
  TYPE :: real_array
    CHARACTER(len=2) :: atom_name
    INTEGER :: time_step
    REAL(kind=8), DIMENSION(3) :: coord_array
    TYPE (real_array), POINTER :: p
  END TYPE real_array

  ! Derived data type to store water molecules around a center atom
  TYPE :: water_array
    INTEGER :: time_step
    CHARACTER(2) :: center_atom_name
    REAL(kind=8), DIMENSION(9) :: coord_water_array
    REAL(kind=8) :: d_RO
    REAL(kind=8), DIMENSION(2) :: theta 
    TYPE (water_array), POINTER :: p
  END TYPE water_array

  ! Define the derived data type for water
  TYPE :: water
    INTEGER :: molecule_id
    TYPE (atom):: oxygen
    TYPE (atom):: hydrogen1
    TYPE (atom):: hydrogen2
    REAL(kind=8) :: mass
  END TYPE water

  ! Define the derived data type for a water pair
  TYPE :: water_pair
    INTEGER :: water_pair_id
    INTEGER :: head_id ! For example, 2
    INTEGER :: tail_id ! For example, 1
    TYPE (water) :: head ! water molecule
    TYPE (water) :: tail  !water molecule
    LOGICAL :: h_state 
    REAL(kind=8) :: h ! Hydrogen bond population operator 1 or 0 
    REAL(kind=8) :: hb_length 
  END TYPE water_pair

INTERFACE

  SUBROUTINE read_surf_traj(surf_filename,nmo_start,nmo_end,ns,n_grid,n_samples,surf_info)
      ! To read info from the isosurface trajectory file (format: surf_traj_*.xyz)
      IMPLICIT NONE
      INTEGER ::ierror, i_sample, i_grid
      INTEGER, INTENT(IN) :: n_grid ! KEEP
      INTEGER, INTENT(IN) :: n_samples ! n_samples = INT(nmo/ns) KEEP
      INTEGER, INTENT(IN) :: ns ! Get one sample from the trajectory every ns step.
      INTEGER, INTENT(IN) :: nmo_start, nmo_end ! To get the total number of moves
      REAL(kind=8), DIMENSION(2,n_grid,n_samples), INTENT(IN) :: surf_info
      CHARACTER(LEN=4) :: head_char
      character(LEN=200) :: surf_filename
      INTEGER :: y
      INTEGER,PARAMETER :: indx = 20
  END SUBROUTINE read_surf_traj

  REAL(KIND=8) FUNCTION distance2(r1,r2,boxsize)
      ! Date: 2021-3-24
      ! This is the correct version for calculate distance between two particles with PBC considered.
      IMPORT :: rk
      IMPLICIT NONE
      REAL(kind=rk), DIMENSION(3),INTENT(IN) :: r1,r2
      REAL(kind=rk), DIMENSION(3),INTENT(IN) :: boxsize
      REAL(kind=rk) :: dx,dy,dz
      REAL(kind=rk) :: temp
  END FUNCTION distance2
  
  REAL(KIND=8) FUNCTION distance2_(u1,v1,w1,u2,v2,w2,a,b,c)
      IMPLICIT NONE
      integer,parameter :: rk=8  
      REAL(kind=rk),INTENT(INOUT) :: u1,v1,w1,u2,v2,w2,a,b,c
      logical :: A1,A2,A3,B1,B2,B3,C1,C2,C3
  END FUNCTION distance2_

  REAL(kind=8) FUNCTION P2(x)
      IMPLICIT NONE
      ! P2: 2-order Legendre polynomial
      integer,parameter :: rk=8  
      real(kind=rk),INTENT(IN) :: x
  END FUNCTION P2 

  REAL(KIND=8) FUNCTION diff_axis_v1(u1,u2,a)
      IMPLICIT NONE
      integer,parameter :: rk=8  
      logical :: A1,A2,A3
      REAL(kind=rk),INTENT(INOUT) :: u1, u2
      REAL(kind=rk),INTENT(IN) :: a
  END FUNCTION diff_axis_v1

  REAL(kind=8) FUNCTION diff_axis(u1,u2,h)
      IMPLICIT NONE
      integer,parameter :: rk=8  
      REAL(kind=rk),INTENT(IN) :: u1, u2
      REAL(kind=rk),INTENT(IN) :: h
      REAL(kind=rk) :: du 
  END FUNCTION diff_axis
  
  REAL(kind=8) FUNCTION pm_adh(r1,r2,r3,boxsize)
      IMPLICIT NONE
      integer,parameter :: rk=8  
      REAL(kind=rk),DIMENSION(3),INTENT(IN) :: r1,r2,r3
      REAL(kind=rk),DIMENSION(3),INTENT(IN) :: boxsize 
      REAL(kind=rk) :: diff_axis
  END FUNCTION pm_adh

  REAL(kind=8) FUNCTION pm_ahd(r1,r2,r3,boxsize)
      IMPLICIT NONE
      integer,parameter :: rk=8  
      REAL(kind=rk),DIMENSION(3),INTENT(IN) :: r1,r2,r3
      REAL(kind=rk),DIMENSION(3),INTENT(IN) :: boxsize 
      REAL(kind=rk) :: diff_axis
  END FUNCTION pm_ahd 
  
  INTEGER FUNCTION get_total_number_of_lines(file_)
      IMPLICIT NONE
      INTEGER :: nlines, io
      CHARACTER (len=*),INTENT(IN) :: file_
  END FUNCTION 

  FUNCTION hydrogen_ndx_list(ndx_oxygen_1, & 
             ndx_oxygen_2, pos_filename, natoms,boxsize) 
      IMPLICIT NONE
      integer, parameter :: rk=8  
      INTEGER, INTENT(IN) :: ndx_oxygen_1, ndx_oxygen_2
      character(LEN=200), INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      REAL(kind=rk), dimension(3), INTENT(IN) :: boxsize
      !Local 
      REAL(kind=rk), dimension(3) :: r_a, r_b
      real(kind=rk), parameter :: r_ohc=1.21d0   ! rOH (1.1**2)
      real(kind=rk) :: r ! distance between O and H.
      integer :: i,j,nmovie,iatom,& 
                    m1,m2,m3,i_H,&
                    i1,i2,ii,jj,i_OW, num
      INTEGER, DIMENSION(4) :: hydrogen_ndx_list
      real(kind=rk), allocatable, dimension(:,:) :: x,y,z
      character(LEN=3), allocatable, dimension(:) :: atom_type
      integer, allocatable, dimension(:) :: ndx_OW,ndx_H
  END FUNCTION hydrogen_ndx_list

  FUNCTION hydrogen_ndx_list_XO(ndx_X, & 
           ndx_oxygen_2,pos_filename, natoms,boxsize) 
      IMPLICIT NONE
      ! In this function, we use ndx_X replace ndx_oxygen_1. 
      ! ndx_X is not bonded with any Hydrogen atoms
      integer,parameter :: rk=8              

      INTEGER,INTENT(IN) :: ndx_X, ndx_oxygen_2
      character(LEN=200),INTENT(IN) :: pos_filename      ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      REAL(kind=rk),dimension(3),INTENT(IN) :: boxsize
      !Local 
      REAL(kind=rk),dimension(3) :: r_a, r_b
      real(KIND=rk),parameter :: r_ohc=1.21d0   ! rOH (1.1**2)
      real(kind=rk) :: r ! distance between O and H.
      integer :: i,j,nmovie,iatom,& 
              m1,m2,m3,i_H,&
              i1,i2,ii,jj,i_O,num
      INTEGER,DIMENSION(2) :: hydrogen_ndx_list_XO ! the length is 2, because there are 2 Hydrogen are bonded to ndx_oxygen_2 but no Hydrogen is bonded to ndx_X
      real(kind=rk),allocatable,DIMENSION(:,:) :: x,y,z
      character(LEN=3),allocatable,DIMENSION(:) :: atom_type
      integer,allocatable,DIMENSION(:) :: ndx_O,ndx_H
  END FUNCTION hydrogen_ndx_list_XO

  FUNCTION hydrogen_ndx_list_for_oxygen(ndx_oxygen_1, & 
           pos_filename,natoms,boxsize) 
      implicit none
      integer,parameter :: rk=8  
      INTEGER,INTENT(IN) :: ndx_oxygen_1
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      real(kind=rk),dimension(3),INTENT(IN) :: boxsize
      !Local 
      REAL(kind=rk) :: distance2
      INTEGER :: get_number_of_oxygen,get_number_of_hydrogen
      integer :: i,j,iatom,nmovie,& 
                    m1,m3,i_H,&
                    i1,ii,jj,i_OW, num
      INTEGER, DIMENSION(2) :: hydrogen_ndx_list_for_oxygen
      real(kind=rk), allocatable,dimension (:,:) :: x,y,z
      character(LEN=3), allocatable, dimension (:) :: atom_type
      integer, allocatable, dimension (:) :: ndx_OW,ndx_H
      real(kind=rk), dimension(3) :: r_a, r_b
      real(kind=rk), parameter :: r_ohc=1.21d0   ! rOH (1.1**2)
      real(kind=rk) :: r ! distance between O and H.
  END FUNCTION hydrogen_ndx_list_for_oxygen

  FUNCTION get_ndx_of_oxygen(pos_filename, natoms) 
      implicit none
      !Parameters 
      INTEGER, PARAMETER :: rk=8

      !variables
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER, DIMENSION(natoms) :: get_ndx_of_oxygen 
      !Local 
      integer :: i,nmovie,iatom,imovie
      real(kind=rk), allocatable, dimension (:,:)  :: x,y,z
      character(LEN=3), allocatable, dimension (:) :: atom_type
  END FUNCTION get_ndx_of_oxygen

  FUNCTION get_number_of_oxygen(pos_filename, natoms) 
      IMPLICIT NONE

      integer,parameter :: rk=8              

      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      INTEGER :: get_number_of_oxygen 
      !Local 
      integer :: i,nmovie,iatom,imovie
      real(kind=rk),allocatable,DIMENSION(:,:) :: x,y,z
      character(LEN=3),allocatable,DIMENSION(:) :: atom_type
  END FUNCTION get_number_of_oxygen

  FUNCTION get_number_of_hydrogen(pos_filename, natoms) 
      IMPLICIT NONE
      integer,parameter :: rk=8              
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      INTEGER :: get_number_of_hydrogen
      !Local 
      integer :: i,nmovie,iatom,imovie
      real(kind=rk),allocatable,DIMENSION(:,:) :: x,y,z
      character(LEN=3),allocatable,DIMENSION(:) :: atom_type
  END FUNCTION get_number_of_hydrogen

  FUNCTION get_number_of_iodine(pos_filename, natoms) 
      IMPLICIT NONE
      integer,parameter :: rk=8              
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      INTEGER :: get_number_of_iodine 
      !Local 
      integer :: i,nmovie,iatom,imovie
      real(kind=rk),allocatable,DIMENSION(:,:) :: x,y,z
      character(LEN=3),allocatable,DIMENSION(:) :: atom_type
  END FUNCTION get_number_of_iodine

  SUBROUTINE ghbond_interface(filename,list,nmo,nat,arr)
      !Purpose: To generate the index list of all Oxygen-Oxygen pairs (quasi-bonds) which are once located in
      !         the interface of neat water. There's an inportant idea in this code: We assume that ergodity is not satisfied,
      !         therefore, for a given layer thickness, we list all O-O pairs' indices in one list file, 
      !         and calculate the hydrogen bond dynamics for all these pairs (quasi-HBs), the correlation of some
      !         pairs may calculated serveral times. It is exactly what we want, because the more the pair is
      !         considered, the more weights this pair will have in the statistis. After all the correlations
      !         are calculated, we calculate the average correlation over all these correlations.
      !Output: A list file of O-O pairs. Note that one O-O pair may occurs more than once, because the O-O
      !        pair may in interface at both time t_i and t_j.
      IMPLICIT NONE
      character(LEN=200),INTENT(IN) :: filename            ! specific filename to analyz data
      character(LEN=200),INTENT(INOUT) :: list
      INTEGER,INTENT(IN) :: nmo ! steps of trajectory
      INTEGER,INTENT(IN) :: nat ! number of atoms
      !To save the indices of the molecules for generating list file, we define an array for each time point (jj, in this code)
      INTEGER,DIMENSION(nmo,nat),INTENT(INOUT) :: arr
      ! Local variables
      INTEGER,PARAMETER :: step=100 ! The parameter 'step' should not be too small, otherwise, you will waste your time on many repeated calculation. 
                                    ! Here, 100 means every 100* ns steps, we select the molecules in interface.
      integer :: i,n,m1,m2,i1,i2,jj
      integer,DIMENSION(nat) :: ndx_O
  END SUBROUTINE ghbond_interface 

  SUBROUTINE ghbond_lino3(filename,pos_filename,natoms,&
      list_ow_ow_pairs,list_on_ow_pairs,boxsize)

      IMPORT distance2
      !Purpose: To generate the index list of possible Nitrate O -- water O pairs, and water O --water O pairs (quasi-bonds)
      !Returns: two list file. One is the Nitrate O -- water O pairs
      !         The other is water O --water O pairs
      !These Nitrete O -- WaterO pairs help one to define Nitrate O -- OW Hydrogen bonds; 
      ! and these Water O -- Water O pairs helps to define OW-OW hydrogen bonds.
      ! This code modified from '/home/gang/Github/hbacf/__hbacf_continuous/__hbond_all_pair_lino3/ghbond_lino3.f95'.
      implicit none
      !==========
      !Parameters
      !==========
      integer,parameter :: rk=8              
      !=========
      !variables
      !=========
      character(LEN=200),INTENT(IN) :: filename            ! specific filename to analyz data
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      character(LEN=200),INTENT(INOUT) :: list_ow_ow_pairs, &
          list_on_ow_pairs
      INTEGER, INTENT(IN) :: natoms
      REAL(kind=rk),dimension(3), INTENT(INOUT) :: boxsize  ! box size
      ! Local variables
      real(kind=rk),parameter :: r_ohc=1.44d0      ! rOH (1.2**2)
      real(kind=rk),parameter :: r_onc=2.25d0      ! rOH (1.5**2)
      real(kind=rk) :: r23
      integer :: i,ii,i4,i5,i6,j,nmovie,iatom,& 
                    imovie,m1,m2,m3,&
                    i_O,i_H,i_N,i_ON,i_OW
      real(kind=rk), allocatable,dimension (:,:) :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:) :: ndx_O,ndx_H,&
          ndx_N,ndx_OW,ndx_ON
  END SUBROUTINE ghbond_lino3 

  INTEGER FUNCTION grid_index(x,y,divx,divy,nb_divy)
      ! transfer the coordinates (x,y) to grid_index, which is an integer
      IMPLICIT NONE
      !==========
      !Parameters
      !==========
      integer,parameter :: rk=8              
      !=========
      !variables
      !=========
      INTEGER,DIMENSION(2) :: ind
      REAL(kind=rk),INTENT(IN) :: divx, divy
      INTEGER,INTENT(IN) :: nb_divy
      REAL(kind=rk),INTENT(IN) :: x,y
  END FUNCTION grid_index

  LOGICAL FUNCTION oh_in_surf1(surf1_mol1,z1,thickness)
      implicit none
      integer, parameter :: rk=8
      LOGICAL :: mol1_in_surf1
      REAL(kind=rk), INTENT(IN) :: surf1_mol1
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z1
  END FUNCTION oh_in_surf1    

  LOGICAL FUNCTION oh_in_surf2(surf2_mol1,z1,thickness)
      implicit none
      integer, parameter :: rk=8
      LOGICAL :: mol1_in_surf2
      REAL(kind=rk), INTENT(IN) :: surf2_mol1
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z1
  END FUNCTION oh_in_surf2    

  LOGICAL FUNCTION pair_in_surf1(surf1_mol1,z1,surf1_mol2,z2,thickness)
      IMPLICIT NONE
      integer,parameter :: rk=8              
      LOGICAL :: mol1_in_surf1, mol2_in_surf1
      REAL(kind=rk),INTENT(IN) :: surf1_mol1, surf1_mol2
      REAL(kind=rk),INTENT(IN) :: thickness
      REAL(kind=rk),INTENT(IN) :: z1, z2
  END FUNCTION pair_in_surf1    

  LOGICAL FUNCTION pair_in_surf2(surf2_mol1,z1,surf2_mol2,z2,thickness)
      IMPLICIT NONE
      integer, parameter :: rk=8              
      LOGICAL :: mol1_in_surf2,mol2_in_surf2
      REAL(kind=8), INTENT(IN) :: surf2_mol1,surf2_mol2
      REAL(kind=8), INTENT(IN) :: thickness
      REAL(kind=8), INTENT(IN) :: z1,z2
  END FUNCTION pair_in_surf2    

  LOGICAL FUNCTION mol_in_surf1(surf1_mol,z1,thickness)
      IMPLICIT NONE
      !Check if a molecule is in the surf1 (the lower layer)
      integer, parameter :: rk=8
      REAL(kind=rk),INTENT(IN) :: surf1_mol
      REAL(kind=rk),INTENT(IN) :: thickness
      REAL(kind=rk),INTENT(IN) :: z1
  END FUNCTION mol_in_surf1    

  LOGICAL FUNCTION mol_in_surf2(surf2_mol,z2,thickness)
      IMPLICIT NONE
      integer,parameter :: rk=8              
      !Check if a molecule is in the surf2 (the upper layer)
      REAL(kind=rk),INTENT(IN) :: surf2_mol
      REAL(kind=rk),INTENT(IN) :: thickness
      REAL(kind=rk),INTENT(IN) :: z2
  END FUNCTION mol_in_surf2    
  
  character(len=20) function int2str(k)
      !  "Convert an integer to string."
      integer, intent(in) :: k
  end function int2str

  character(len=20) function str(k)
      !  "Convert an integer/real to string."
      real (kind=8), intent(in) :: k
  end function str
  
  FUNCTION nth(k,n)
      !To return a string containing the first N characters
      !of the alphabet.
      !Declaring calling parameters:
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: n            ! Length of string to return
      CHARACTER(len=20),INTENT(IN) :: k  ! String which is adjustl-ed
      CHARACTER(len=n) nth                ! Returned string
  END FUNCTION nth

  INTEGER FUNCTION sampling_number(nmo_start,nmo_end,ns)
      !To calculate the total numbers of samples one want to include in their analysis.
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ns ! Get one sample from the trajectory every ns step.
      INTEGER,INTENT(IN) :: nmo_start, nmo_end ! To get the total number of moves
  END FUNCTION sampling_number

  SUBROUTINE read_traj(indx,nmo_start,nmo_end,ns,nat,n_samples,sampled_movie,sampled_time,atom_info)
      ! To read info from the trajectory file (format: ***.xyz)
      ! to READ data starting from a pattern-matched line.
      IMPORT :: atom
      IMPLICIT NONE
      integer,parameter :: rk=8              
      CHARACTER(LEN=4) :: head_char
      INTEGER :: iatom,i_sample
      INTEGER,INTENT(IN) :: indx
      INTEGER,INTENT(IN) :: nat
      INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
      INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
      INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.

      TYPE(atom),DIMENSION(nat,n_samples),INTENT(INOUT) :: atom_info
      INTEGER,DIMENSION(n_samples) :: sampled_movie
      REAL(kind=rk),DIMENSION(n_samples) :: sampled_time
      INTEGER :: y
  END SUBROUTINE read_traj

  SUBROUTINE read_traj_sphere(indx,nmo_start,nmo_end,ns,nat,n_samples,sampled_movie,sampled_time,sampled_energy,sphere_info)
      ! To read info from the trajectory file (format: ***.xyz)
      ! to READ data starting from a pattern-matched line.
      IMPORT :: sphere
      IMPLICIT NONE
      integer,parameter :: rk=8              
      CHARACTER(LEN=4) :: head_char
      INTEGER :: iatom,i_sample
      INTEGER,INTENT(IN) :: indx
      INTEGER,INTENT(IN) :: nat
      INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
      INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
      INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.

      TYPE(sphere),DIMENSION(nat,n_samples),INTENT(INOUT) :: sphere_info
      INTEGER,DIMENSION(n_samples) :: sampled_movie
      REAL(kind=rk),DIMENSION(n_samples) :: sampled_time, sampled_energy
      INTEGER :: y
  END SUBROUTINE read_traj_sphere

  SUBROUTINE read_interface_input(boxsize,delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,&
               criterion,surf_filename,thickness) ! RIGHT one 
      IMPLICIT NONE
      integer, parameter :: rk=8              
      REAL(kind=rk), dimension(3), INTENT(INOUT) :: boxsize  
      integer, INTENT(INOUT) :: criterion
      REAL(kind=rk), INTENT(INOUT) :: delta_t0  ! For reading data
      character(LEN=200), INTENT(INOUT) :: filename
      INTEGER :: i,iatom, ierr, imovie
      character(LEN=50) :: line 
      INTEGER, INTENT(INOUT) :: nat ! number of atoms
      INTEGER, INTENT(INOUT) :: nmo_start 
      INTEGER, INTENT(INOUT) :: nmo_end
      INTEGER, INTENT(INOUT) :: ns
      character(LEN=200), INTENT(INOUT) :: pos_filename
      character(LEN=200), INTENT(INOUT) :: surf_filename
      REAL(KIND=rk), INTENT(INOUT) :: thickness ! the thickness of the instantaneous interfaces
  END SUBROUTINE read_interface_input

  SUBROUTINE read_surf_coord(indx,n_samples,n_grid,surf_info_fortran)
    IMPLICIT NONE
    INTEGER :: ierror, i_sample, i_grid
    INTEGER :: skip_num ! The number of steps skipped before reading
    INTEGER,INTENT(IN) :: indx
    INTEGER,INTENT(IN) :: n_grid ! n_grid = nb_divx * nb_divy; KEEP
    INTEGER,INTENT(IN) :: n_samples ! n_samples = INT(nmo/ns); KEEP
    REAL(kind=8),DIMENSION(2,n_samples,n_grid),INTENT(INOUT) :: surf_info_fortran
  END SUBROUTINE read_surf_coord

  SUBROUTINE load(pos_filename,nmo_start,nmo_end,nat,ns,n_samples,boxsize, &
           sampled_pos_filename,sampled_movie,sampled_time, &
           nb_divx,nb_divy,nb_divz,n_grid,divx,divy,divz,whish_size,atom_info)
      IMPORT :: atom
      IMPLICIT NONE
      integer, parameter :: rk=8

      REAL(kind=rk), dimension(3),INTENT(INOUT) :: boxsize
      REAL(kind=rk), INTENT(INOUT) :: divx, divy, divz
      INTEGER, INTENT(IN) :: nat ! number of atoms or nb_atoms
      INTEGER, INTENT(INOUT) :: nb_divx, nb_divy, nb_divz, n_grid 
      INTEGER, INTENT(IN) :: nmo_end  
      INTEGER, INTENT(IN) :: nmo_start  
      INTEGER, INTENT(IN) :: ns ! Get one sample from the trajectory every ns step.
      INTEGER, INTENT(IN) :: n_samples ! n_samples = INT(nmo/ns)
      character(LEN=*),INTENT(IN) :: pos_filename

      INTEGER, DIMENSION(n_samples),INTENT(INOUT) :: sampled_movie
      TYPE(atom), DIMENSION(nat,n_samples),INTENT(INOUT) :: atom_info !Should be declared after n_samples is decared
      character(LEN=*), INTENT(INOUT):: sampled_pos_filename
      REAL(kind=rk), DIMENSION(n_samples), INTENT(INOUT) :: sampled_time
      REAL(kind=rk), INTENT(IN) :: whish_size ! Angstrom

      !Local varables
      INTEGER,PARAMETER :: n_axis = 2 ! n_axis=2 means z-axis is the normal axis.
      REAL(kind=rk) :: center_pos 
      INTEGER :: i,iatom,imovie
      INTEGER :: num_wat_pairs 
      REAL(kind=rk) :: sum_mass
  END SUBROUTINE load 

  SUBROUTINE molecules_in_interface(n_samples,nat,arr,atom_info, &
         n_grid,divx,divy,divz,nb_divx,nb_divy,nb_divz,thickness,&
         surf_info)
      !=========================================================
      !Purpose: to obtain the indices of molecules in interfaces
      !         based on selection of molecules. That is, only in
      ! Sulpizi's method, we need this function 
      ! 'molecules_in_interface()', in which surf_info is an 
      ! array.
      ! Output: An 2D array 'indx_array': 
      !         the first axis: time labels (jj)
      !         the 2nd axis: atom's labels (m) 
      !=========================================================
      IMPORT :: atom
      IMPLICIT NONE

      integer, parameter :: rk=8 ! local 

      INTEGER,INTENT(IN) :: nat ! number of atoms
      INTEGER,INTENT(IN) :: n_samples ! steps of trajectory, ie, n_samples 
      INTEGER,INTENT(IN) :: nb_divx,nb_divy,nb_divz,n_grid 
      REAL(kind=rk),INTENT(IN) :: divx,divy,divz
      REAL(KIND=rk), INTENT(IN) :: thickness
      !Local variables
      LOGICAL :: condition1,condition2  
      INTEGER :: index_mol
      INTEGER :: jj,m,n 
      !To save the indices of the molecules for generating list file, we define a array for each time point (jj, in this code)
      INTEGER,DIMENSION(n_samples,nat),INTENT(INOUT) :: arr
      TYPE(atom),DIMENSION(nat,n_samples),INTENT(IN) :: atom_info
      REAL(kind=rk), DIMENSION(2,n_grid,n_samples), INTENT(IN) :: surf_info
  END SUBROUTINE molecules_in_interface 

  SUBROUTINE ghbacf_interface_c_pbc_format2(boxsize,delta_t0,&
          filename,pos_filename,list_filename,n_samples,nat,ns,&
          criterion,atom_info,n_grid,divx,divy,divz,nb_divx,nb_divy,&
          nb_divz,thickness,surf_info)
      IMPORT :: atom
      IMPLICIT NONE
      !I HAS USED: get_total_numer_of_lines,&
      !          hydrogen_ndx_list,&
      !          distance2,pm_adh,pm_ahd,& 
      !          grid_index,pair_in_surf1,pair_in_surf2,&
      !          str,nth
      integer,parameter :: rk=8 ! local 

      character(LEN=200),INTENT(INOUT) :: filename,pos_filename
      character(LEN=200),INTENT(IN) :: list_filename
      integer,INTENT(IN) :: criterion
      INTEGER,INTENT(IN) :: nat ! number of atoms
      INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
      REAL(kind=rk),dimension(3),INTENT(IN) :: boxsize
      TYPE(atom),DIMENSION(nat,n_samples),INTENT(IN) :: atom_info
      REAL(kind=rk),INTENT(IN) :: thickness ! the thickness of the instantaneous interfaces
      INTEGER,INTENT(IN) :: nb_divx, nb_divy, nb_divz, n_grid 
      REAL(kind=rk),INTENT(IN) :: divx, divy, divz
      REAL(KIND=rk),DIMENSION(2,n_grid,n_samples),INTENT(IN) :: surf_info
      
      !Local variables
      CHARACTER(len=1) :: char_thickness ! for saving the thickness in the files' names
      LOGICAL :: condition1, condition2
      real(kind=rk),parameter :: cosPhiC123=0.866d0 ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: cosPhiC132=-0.5d0 ! -1./2; phiC132=2pi/3.
      REAL(kind=rk),parameter :: h_min=0.5d0 ! condition for the existence of a h-bond for a step
      REAL(kind=rk),parameter :: hb_min=0.5d0 ! condition for the existence of h-bond for a pair of water molecules
      REAL(kind=rk) :: r13,cosphi,pm, cosphi_, pm_
      REAL(kind=rk) :: r21,r31,r32,r23 ! For the second criterion of HB
      REAL(kind=rk) :: qj,tot_hb,delta_t,delta_t0,hb_per_frame,ave_h
      REAL(kind=rk),dimension(3) :: r1, r2, r3 ! pbc 
      integer :: m1,m2,m3,mt,nqj,tot_nhb,n_bonded_pairs,ns
      REAL(kind=rk),allocatable,DIMENSION(:) :: h,hb,corr_h
      real(kind=rk),allocatable,DIMENSION(:,:) :: x,y,z
      integer,allocatable,dimension(:) :: ndx_1,ndx_2,nhb_exist
      integer,dimension(4) :: ndx_3_list
      REAL(kind=rk) :: scalar 
      logical,allocatable,DIMENSION(:) :: hb_exist
      INTEGER :: nmo  ! nmo is not necessary, we set nmo = n_samples, because we do not want to change too much
      INTEGER :: nwat ! number of water molecules
      INTEGER :: i,j,jj,k 
      INTEGER :: index_mol1, index_mol2
      real(kind=rk),parameter :: rooc=12.25d0 ! cutoff distance of rOO (3.5**2 )
  END SUBROUTINE 

  SUBROUTINE ghbacf_interface_n_pbc_format2(boxsize,delta_t0, &
      filename,pos_filename,list_filename,n_samples,nat,ns, &
      criterion,atom_info,n_grid,divx,divy,divz,nb_divx,nb_divy,&
      nb_divz,thickness,surf_info)
      IMPORT :: atom
      IMPLICIT NONE
      character(LEN=200),INTENT(INOUT) :: filename,pos_filename
      character(LEN=200),INTENT(IN) :: list_filename
      integer,INTENT(IN) :: criterion
      INTEGER,INTENT(IN) :: nat ! number of atoms
      INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
      integer,parameter :: rk=8 ! local 
      REAL(kind=rk),dimension(3),INTENT(IN) :: boxsize
      TYPE(atom),DIMENSION(nat,n_samples),INTENT(IN) :: atom_info
      REAL(kind=rk),INTENT(IN) :: thickness ! the thickness of the instantaneous interfaces
      INTEGER,INTENT(IN) :: nb_divx,nb_divy,nb_divz,n_grid 
      REAL(kind=rk),INTENT(IN) :: divx, divy,divz
      REAL(KIND=rk),DIMENSION(2,n_grid,n_samples),INTENT(IN) :: surf_info
      !Local variables
      real(kind=rk),parameter :: rooc=12.25d0                 ! cutoff distance of rOO (3.5**2 )
      real(kind=rk),parameter :: cosPhiC123=0.866d0           ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: cosPhiC132=-0.5d0            ! -1./2; phiC132=2pi/3.
      REAL(kind=rk),parameter :: h_min=0.5d0  ! condition for the existence of a h-bond for a step
      REAL(kind=rk),parameter :: hb_min=0.5d0 ! condition for the existence of h-bond for a pair of water molecules
      REAL(kind=rk) :: r13, cosphi, pm, cosphi_, pm_
      REAL(kind=rk) :: r21, r31, r32, r23 ! For the second criterion of HB
      REAL(kind=rk) ::qj, tot_hb, delta_t, delta_t0, hb_per_frame, ave_h
      REAL(kind=rk),dimension(3) :: r1, r2, r3 ! pbc 
      integer :: m1,m2,m3,mt,nqj,tot_nhb,n_bonded_pairs,ns
      REAL(kind=rk),allocatable,DIMENSION(:)  :: h,h_d,hb,corr_n
      real(kind=rk),allocatable,DIMENSION(:,:) :: x,y,z
      integer,allocatable,dimension(:) :: ndx_1,ndx_2,nhb_exist
      integer,dimension(4) :: ndx_3_list
      REAL(kind=rk) :: scalar 
      logical,allocatable,DIMENSION(:)  :: hb_exist
      INTEGER :: nmo  ! nmo is not necessary, we set nmo = n_samples, because we do not want to change too much
      INTEGER :: nwat ! number of water molecules
      INTEGER :: i,j,k,jj 
      CHARACTER(len=1) :: char_thickness ! for saving the thickness in the files' names
      INTEGER :: index_mol1,index_mol2
      LOGICAL :: condition1,condition2
  END SUBROUTINE 

  SUBROUTINE ghbacf_interface_k_pbc_format2(boxsize,delta_t0, & 
      filename,pos_filename,list_filename,n_samples,nat, &
      ns,criterion,atom_info,n_grid,divx,divy,divz,nb_divx, &
      nb_divy,nb_divz,thickness,surf_info)
      IMPORT :: atom
      IMPLICIT NONE
      character(LEN=200),INTENT(INOUT) :: filename,pos_filename
      character(LEN=200),INTENT(IN) :: list_filename
      integer,INTENT(IN) :: criterion
      INTEGER,INTENT(IN) :: nat ! number of atoms
      INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
      integer,parameter :: rk=8 ! local 
      REAL(kind=rk),dimension(3),INTENT(IN) :: boxsize
      REAL(kind=rk),INTENT(IN) :: delta_t0
      TYPE(atom),DIMENSION(nat,n_samples),INTENT(IN) :: atom_info
      REAL(kind=rk),INTENT(IN) :: thickness ! the thickness of the instantaneous interfaces
      INTEGER,INTENT(IN) :: nb_divx, nb_divy, nb_divz, n_grid 
      REAL(kind=rk),INTENT(IN) :: divx, divy, divz
      REAL(KIND=rk),DIMENSION(2,n_grid,n_samples),INTENT(IN) :: surf_info
      !Local variables
      real(kind=rk),parameter :: rooc=12.25d0                 ! cutoff distance of rOO (3.5**2 )
      real(kind=rk),parameter :: cosPhiC123=0.866d0              ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: cosPhiC132=-0.5d0            ! -1./2; phiC132=2pi/3.
      REAL(kind=rk),parameter :: h_min=0.5d0 ! condition for the existence of a h-bond for a step
      REAL(kind=rk),parameter :: hb_min=0.5d0 ! condition for the existence of h-bond for a pair of water molecules
      REAL(kind=rk) :: r13,cosphi,pm, cosphi_, pm_
      REAL(kind=rk) :: r21,r31,r32,r23 ! For the second criterion of HB
      REAL(kind=rk) :: qj,tot_hb,delta_t,ddelta_t,hb_per_frame,ave_h
      REAL(kind=rk),dimension(3) :: r1,r2,r3 ! pbc 
      integer :: m1,m2,m3,mt,nqj,tot_nhb,n_bonded_pairs,ns
      REAL(kind=rk),allocatable,DIMENSION(:) :: h,hb,corr_h,dc
      real(kind=rk),allocatable,DIMENSION(:,:) :: x,y,z
      integer,allocatable,dimension(:) :: ndx_1,ndx_2,nhb_exist
      integer,dimension(4) :: ndx_3_list
      CHARACTER(len=1) :: char_thickness ! for saving the thickness in the files' names
      REAL(kind=rk) :: scalar 
      logical,allocatable,DIMENSION(:) :: hb_exist
      INTEGER :: nmo ! nmo is not necessary, we set nmo = n_samples, because we do not want to change too much
      INTEGER :: nwat ! number of water molecules
      INTEGER :: i,j,k,jj 
      INTEGER :: index_mol1,index_mol2
      LOGICAL :: condition1,condition2
  END SUBROUTINE 

  SUBROUTINE skip_lines(indx,i_input)
      ! To skip lines when read data from the input
      IMPLICIT NONE
      INTEGER :: i
      INTEGER,INTENT(IN) :: i_input,indx
  END SUBROUTINE skip_lines

  SUBROUTINE make_neighbors(nat,n_samples,boxsize,neighbors_number,neighbors_indices,sphere_info)
    IMPORT :: distance2,sphere,num_coord_max,rk,r_OH_LU2,r_ON_LU2
    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER,INTENT(IN) :: nat, n_samples
    real(kind=rk), DIMENSION(3), INTENT(IN) :: boxsize
    TYPE(sphere),DIMENSION(nat,n_samples),INTENT(INOUT) :: sphere_info
    INTEGER,DIMENSION(nat),INTENT(INOUT) :: neighbors_number
    INTEGER,DIMENSION(nat,num_coord_max),INTENT(INOUT) :: neighbors_indices
    !REAL(kind=rk) :: distance2

    INTEGER,PARAMETER :: start_step = 1
    REAL(kind=rk) :: r_LU2  ! LU: Least Upper bound; 2: squared value
  END SUBROUTINE make_neighbors

  SUBROUTINE obtain_indices(atom_name,indices,nat,n_samples,sphere_info)
    !===================================================================
    ! In this subroutine, I want to do one thing for an chemical element (X, a string denoted as atom_name):
    ! obtain the indices of atom X and store it in the array indices
    ! indices: the output array
    !===========================
    IMPORT :: sphere, num_coord_max
    IMPLICIT NONE

    CHARACTER(len=*) :: atom_name 
    INTEGER,INTENT(IN) :: nat,n_samples
    TYPE(sphere),DIMENSION(nat,n_samples),INTENT(INOUT) :: sphere_info
    INTEGER,DIMENSION(nat), INTENT(INOUT) :: indices
    
    INTEGER,PARAMETER :: start_step=1 
    INTEGER :: i, n_count
  END SUBROUTINE obtain_indices

  SUBROUTINE make_wat_neighbors_for_X_shell(atom_name,step_t,nat,n_samples,boxsize, &
             wat_neighbors_number,wat_neighbors_indices,sphere_info)
    !===================================================================
    ! In this subroutine, I want to do one thing:
    ! For an atom (index = i_ion), at time t (denoted by step_t), 
    ! 1: to find the number (N) of its neighboring water molecules and store N in wat_neighbors_number(i_ion,t); 
    ! 2: to add the index of O of this water molecule into the wat_neighbors_indices(i_ion,t)
    !====================================================================
    
    IMPORT :: rk,sphere,num_coord_mol_max,r_OH_LU2,r_NN_OW_LU2,r_ON_OW_LU2,r_Li_OW_LU2,r_Na_OW_LU2,r_K_OW_LU2,&
              r_OW_OW_LU2,r_I_OW_LU2,distance2
    IMPLICIT NONE

    INTEGER :: i, j, n_count
    CHARACTER(len=*),INTENT(IN) :: atom_name 
    INTEGER,INTENT(IN) :: nat,n_samples
    real(kind=rk), DIMENSION(3), INTENT(IN) :: boxsize
    TYPE(sphere),DIMENSION(nat,n_samples),INTENT(INOUT) :: sphere_info
    INTEGER,DIMENSION(nat),INTENT(INOUT) :: wat_neighbors_number
    INTEGER,DIMENSION(nat,num_coord_mol_max),INTENT(INOUT) :: wat_neighbors_indices ! The max coord. number of molecules 
     
    ! Local variable
    INTEGER,DIMENSION(nat) :: indices
    !REAL(kind=rk) :: dist2
    
    INTEGER, INTENT(IN):: step_t 
    REAL(kind=rk) :: r_LU2 ! LU: Least Upper bound; 2: squared value 
  END SUBROUTINE make_wat_neighbors_for_X_shell

  CHARACTER(LEN=200) FUNCTION name_list_OH_at_shell_interface(filename,core_atom_str)
    IMPLICIT NONE

    CHARACTER(len=*),INTENT(IN) :: filename 
    character(LEN=*),INTENT(IN) :: core_atom_str        ! For example, core_atom_str can be "I","Na", etc.
  END FUNCTION

  CHARACTER(LEN=200) FUNCTION name_c2_at_shell_interface(filename,core_atom_str)
    IMPLICIT NONE

    CHARACTER(len=*),INTENT(IN) :: filename 
    character(LEN=*),INTENT(IN) :: core_atom_str        ! For example, core_atom_str can be "I","Na", etc.
  END FUNCTION
END INTERFACE

END MODULE module_ihb
