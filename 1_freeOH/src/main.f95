PROGRAM main_freeoh
  ! 0) Purpose:     
  !    Determine the free OH density dynamics of interface layer of neat water, by using Sulpizi's selection algorithm.
  !
  ! 0a) Date     Programmer     version
  !  ======      ==========     =======
  !  24/9/8      Gang Huang     Original code
  !============
  ! Declaration
  !============
  USE module_ihb
  USE atom_module, ONLY: atom, atom_info
  USE water_molecule_types
  IMPLICIT NONE

  !==========
  !parameters
  !==========             
  INTEGER, PARAMETER :: nmo_sampling_start= 100 ! Starting step index for starting choosing O-O pairs (This parameter is ONLY needed for case 1.)

  ! The array atom_info can be shared by subroutines  
  !To declear data for calculating the total time cosummed.
  INTEGER :: begin_time,end_time,rat
  !To declare data to share between routines.
  CHARACTER(LEN=200) :: sampled_pos_filename ! Jie: not used
  INTEGER, ALLOCATABLE, DIMENSION(:) :: sampled_movie
  REAL(KIND=rk), ALLOCATABLE, DIMENSION(:) :: sampled_time
  ! For griding
  INTEGER :: nb_divx, nb_divy, nb_divz, n_grid, i_grid
  REAL(KIND=rk) :: divx, divy, divz

  REAL(KIND=rk) :: thickness ! the thickness of the instantaneous interfaces
  REAL(KIND=rk), ALLOCATABLE, DIMENSION(:,:,:) :: surf_info 

  CHARACTER(LEN=2) :: atom_type
  REAL(KIND=rk), dimension(3) :: boxsize
  INTEGER :: criterion
  REAL(kind=rk) :: delta_t0 ! For reading data
  CHARACTER(LEN=200) :: filename
  CHARACTER(LEN=2) :: guest_atom
  CHARACTER(LEN=2) :: host_atom
  INTEGER :: i, iatom, imovie, i_sample
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: indx_array1, indx_array2
  CHARACTER(LEN=200) :: list_oxygen_pairs
  INTEGER :: nat ! number of atoms
  INTEGER :: nmo_start ! starting step index
  INTEGER :: nmo_end ! end step index
  INTEGER :: ns ! Get one sample from the trajectory every ns step. 
  !(NOTE: The ns should be the same as the ns in the 1_chandler_fast.py)
  INTEGER :: ns_2nd ! Jie: 1
  INTEGER :: n_samples ! n_samples = INT(nmo/ns)
  INTEGER :: nwat ! number of water molecules
  CHARACTER(LEN=200) :: surf_filename
  CHARACTER(LEN=200) :: pos_filename
  INTEGER :: n_H ! number of OH groups; or num of H atoms
  INTEGER :: n_O ! number of O atoms
 
  !==============
  !Initialization
  !==============
  iatom=0;imovie=0;i=0
  atom_type=''
  guest_atom="H"
  host_atom="O"
  boxsize=(/0.d0,0.d0,0.d0/)
  nat=0
  nwat=nat/3
  criterion=1 ! 1 means ADH criterion of H-Bbond definition
  filename=""; pos_filename=""; sampled_pos_filename=""
  surf_filename=""
  list_oxygen_pairs=""
  thickness=1.d0
  ns_2nd = 1

  CALL system_clock(begin_time,rat) !To get the starting time
  


  !============================================
  ! To read the required controlling parameters
  !============================================
  CALL read_interface_input(boxsize,delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,&
          criterion,surf_filename,thickness) 

  write(*,*) "boxsize= ", boxsize(1), boxsize(2), boxsize(3)
  nb_divx = nint(boxsize(1)/whish_size) ! round the argument to the nearest integer.
  nb_divy = nint(boxsize(2)/whish_size) ! round the argument to the nearest integer.
  nb_divz = nint(boxsize(3)/whish_size) ! round the argument to the nearest integer.
  divx = boxsize(1)/REAL(nb_divx,rk)
  divy = boxsize(2)/REAL(nb_divy,rk)
  divz = boxsize(3)/REAL(nb_divz,rk)
  n_grid = nb_divx * nb_divy

  write(*,*) "filename= ", filename
  n_samples = sampling_number(nmo_start, nmo_end,ns)
  ALLOCATE(sampled_movie(n_samples))
  ALLOCATE(sampled_time(n_samples))
  n_H = nat * 2/3
  n_O = nat / 3
  ALLOCATE(atom_info(nat,n_samples))
  ALLOCATE(H_info(n_H,n_samples))
  ALLOCATE(O_info(n_O,n_samples))


  !=======================
  !read in trajectory file
  !=======================
  OPEN(10,file=trim(pos_filename))
  CALL read_traj_3(10,nmo_start,nmo_end,ns,nat,n_samples,sampled_movie,sampled_time,boxsize, criterion)
  CLOSE(10)
  WRITE(6,*) 'End of trajectory reading.'

  !====================
  !read surf trajectory
  !====================
  ALLOCATE(surf_info(2,n_grid,n_samples))
  surf_info = 0
  ns_2nd = 1 ! sample freq is 1, ie., all data are sampled
  CALL read_surf_traj(surf_filename,nmo_start,nmo_end,ns_2nd,n_grid,n_samples,surf_info)
  ! Use array instead of linked list, it may be faster. 
  ALLOCATE(indx_array1(n_samples,nat))
  ALLOCATE(indx_array2(n_samples,nat))
  indx_array1 = 0
  indx_array2 = 0
  CALL molecules_in_interface(n_samples,nat,indx_array1,indx_array2,atom_info,&
     n_grid,divx,divy,divz,nb_divx,nb_divy,nb_divz,thickness,surf_info)

  CALL ghbond_interface(filename,list_oxygen_pairs,n_samples,nat,indx_array1,indx_array2)
  !Calculate nf(t)
  CALL ghbacf_interface_nf_pbc_format2(boxsize,delta_t0,filename,pos_filename, &
                    n_samples,nat,ns,criterion,atom_info,n_grid,divx,divy,divz,nb_divx,nb_divy,&
                    nb_divz,thickness,surf_info,O_info,H_info)
  call system_clock(end_time,rat)
  WRITE(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 

END PROGRAM main_freeoh
