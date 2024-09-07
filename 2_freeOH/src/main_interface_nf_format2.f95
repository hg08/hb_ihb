PROGRAM main_interface_nf 
  ! 0) Purpose:     
  ! To calculate the Correlation function C_nf(t) = <nf(0)nf(t)>/<nf>, etc for pure water interface.
  ! If a OH group is free at time t, and is located in one of interfaces, then nf(t)=1, else nf(t)=0
  ! 0a) Date     Programmer     version
  !  ======      ==========     =======
  !  2024/8/31   Gang Huang     Original code
  !============
  ! Declaration
  !============
  USE parameter_shared, ONLY: n_grid
  USE atom_module
  USE water_molecule_types
  USE count_time, ONLY: begin_time,end_time,rat
  USE types
  USE surf_traj
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  integer,parameter :: rk=8 
  INTEGER :: iatom, imovie
  INTEGER :: ns  ! Get one sample from the trajectory every ns step.
  CHARACTER(LEN=2) :: atom_type
  INTEGER :: nmo_start
  INTEGER :: i
  INTEGER :: nmo_end
  INTEGER :: nat  ! number of atoms
  INTEGER :: nwat ! number of water molecules
  INTEGER :: n_samples   ! n_samples = INT(nmo/ns)
  REAL(KIND=rk) :: delta_t0   ! For reading data
  character(LEN=200) :: surf_filename
  character(LEN=200) :: filename, pos_filename
  character(LEN=200) :: list_hydrogen_atoms 
  CHARACTER(LEN=2) :: guest_atom, host_atom
  real(KIND=rk), dimension(3) :: boxsize
  integer :: criterion
  !===============
  ! Initialization
  !===============
  iatom=0;imovie=0;i=0
  atom_type=''
  guest_atom="H"
  host_atom="O"
  boxsize=(/0,0,0/)
  nat=0
  nwat=nat/3
  criterion=1 ! 1 means ADH criterion of H-Bbond definition
  filename=""; pos_filename=""
  surf_filename=""
  !list_oxygen_pairs=""
  list_hydrogen_atoms=""
  call system_clock(begin_time,rat) !To get the starting time
  
  !============================================
  ! To read the required controlling parameters
  !============================================
  CALL read_interface_input(boxsize,delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,&
           criterion,surf_filename) 
  !========================
  ! Sampling the trajectory
  !========================
  CALL load(pos_filename,nmo_start,nmo_end,nat,ns,n_samples,boxsize)
  ! After running the load() subroutine, we have a trajectory file. 
  
  !====================
  !read surf trajectory
  !====================
  open(20,file=trim(surf_filename))
      CALL read_surf_traj(20,nmo_start,nmo_end,ns,n_grid,n_samples)
  close(20)

  CALL gfreeoh(filename,pos_filename,list_hydrogen_atoms,nat)
  CALL ghbacf_interface_nf_pbc_format2(boxsize,delta_t0,filename,&
           n_samples,nat,ns,criterion)
  call system_clock(end_time,rat)
  write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 

END PROGRAM
