PROGRAM main_interface_relax
  ! 0) Purpose:     
  !    Determine the C2(t) dynamics of interface layer of neat water, using Sulpizi's selection algorithm.
  !
  ! Date     Programmer     version
  !======   ==========     =======
  !  21/3/21   Gang Huang     Original code
  !============
  ! Declaration
  !============
  USE module_ihb
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  !The array atom_info can be shared by subroutines  
  TYPE(atom), ALLOCATABLE, DIMENSION(:,:) :: atom_info
  !To declear data for calculating the total time cosummed.
  integer :: begin_time,end_time,rat
  !To declare data to share between routines.
  character(LEN=200) :: sampled_pos_filename
  INTEGER, ALLOCATABLE, DIMENSION(:) :: sampled_movie
  REAL(KIND=rk), ALLOCATABLE, DIMENSION(:) :: sampled_time
  ! For griding
  INTEGER :: nb_divx, nb_divy, nb_divz, n_grid 
  REAL(KIND=rk) :: divx, divy, divz ! divx is how long is each element along x; ...

  REAL(KIND=rk) :: thickness ! the thickness of the instantaneous interfaces
  REAL(KIND=rk),ALLOCATABLE, DIMENSION(:,:,:) :: surf_info_fortran 

  CHARACTER(LEN=2) :: atom_type
  REAL(KIND=rk),dimension(3) :: boxsize
  REAL(kind=rk) :: delta_t0 ! For reading data
  REAL(KIND=rk) :: delta_t  
  character(LEN=200) :: filename
  CHARACTER(LEN=2) :: guest_atom
  CHARACTER(LEN=2) :: host_atom
  INTEGER :: i, iatom, imovie
  !To save the indices of the molecules for generating list file, we define an array for each time point (jj, in this code)
  INTEGER,ALLOCATABLE,DIMENSION(:,:) :: indx_array
  character(LEN=200) :: list_filename
  INTEGER :: nat ! number of atoms
  INTEGER :: nmo_start ! starting step index
  INTEGER :: nmo_end ! end step index
  INTEGER :: ns ! Get one sample from the trajectory every ns step.
  INTEGER :: ns_2nd
  INTEGER :: n_samples ! n_samples = INT(nmo/ns)
  INTEGER :: nwat ! number of water molecules
  character(LEN=200) :: surf_filename
  character(LEN=200) :: pos_filename
  
  !==============
  !Initialization
  !==============
  iatom=0;imovie=0;i=0
  atom_type=" "
  guest_atom="H"
  host_atom="O"
  boxsize=(/0.d0,0.d0,0.d0/)
  nat=0
  nb_divx = 0; nb_divy = 0; nb_divz = 0
  divx = 0; divy = 0; divz = 0
  n_grid = 0
  nwat=nat/3
  filename=" "; pos_filename=" "; sampled_pos_filename=" "
  surf_filename=" "
  list_filename=" "
  thickness=1.d0
  ns_2nd = 1
  delta_t0 = 0.d0
  delta_t = 0.d0  
  call system_clock(begin_time,rat) !To get the starting time
  
  !===========================================
  !To read the required controlling parameters
  !===========================================
  CALL read_interface_input_2(boxsize,delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,&
          surf_filename,thickness) 

  ! Obatin n_samples
  n_samples = sampling_number(nmo_start, nmo_end,ns)
  delta_t = delta_t0 * dble(ns)
  write(*,*) "delta_t = ", delta_t, "(ps)"
  allocate(sampled_movie(n_samples))
  allocate(sampled_time(n_samples))
  allocate(atom_info(nat,n_samples))
 
  !========================
  ! Sampling the trajectory
  !========================
  !!CASE1: If one does not need to recenter, one can just call sample_format2()
  !CALL sample_format2(pos_filename,nmo_start,nmo_end,nat,ns,n_samples)
  !CASE2: If one have to recenter, one call sample_and_recenter_format2() instead.
  CALL sample_and_recenter_format2(pos_filename,nmo_start,nmo_end,nat,ns,n_samples,boxsize,&
       sampled_pos_filename,sampled_movie,sampled_time, &
       nb_divx,nb_divy,nb_divz,n_grid,divx,divy,divz,whish_size,atom_info)
  !WRITE(*,*)"MAIN nb_divx =", nb_divx
  ! After running the sample() or sample_format2() subroutine, therefore we have a new sampled trajectory file (stored in atom_info), 
  ! which is generally shorter than the original one.
  
  !====================
  !read surf trajectory
  !====================
  allocate(surf_info_fortran(2,n_samples,n_grid))
  surf_info_fortran = 0
  !write(*,*) "SHAPE(surf_info_fortran)= ", SHAPE(surf_info_fortran)
  ns_2nd = 1 ! sample freq is 1, ie., all data are sampled
  CALL read_surf_traj(surf_filename,nmo_start,nmo_end,ns,n_grid,n_samples,surf_info_fortran)
  
  ! Use array instead of linked list, it may be faster. 
  allocate(indx_array(n_samples,nat))
  indx_array = 0

  CALL molecules_in_interface(n_samples,nat,indx_array,atom_info,&
     n_grid,divx,divy,divz,nb_divx,nb_divy,nb_divz,thickness,surf_info_fortran)
  !WRITE(*,*)"MAIN4 nb_divx =", nb_divx
  !To determine the indices of Oxygens' pairs that located in one of the interface.
  !CALL ghbond_interface(filename,list_oxygen_pairs,n_samples,nat,indx_array)
  !Similarly, we generate a list for O atoms at one of the interface.
  CALL oh_interface(boxsize,filename,sampled_pos_filename,list_filename,n_samples,nat,indx_array)

  !Calculate n_HB(t),k(t),etc for pure water system. If the format of data is different, one may use another funtion, eg ghbacf_n_pbc_format2().
  CALL relax21(boxsize,filename,sampled_pos_filename,list_filename,delta_t,nat,n_samples,&
      n_grid,divx,divy,divx,nb_divx,nb_divy,nb_divz,thickness,surf_info_fortran) 
  call system_clock(end_time,rat)
  write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 

END PROGRAM main_interface_relax 
