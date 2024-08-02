PROGRAM main_interface_n 
  ! 0) Purpose:     
  ! To calculate the Correlation function C(t) = <h(0)H(t)>/<h>, etc for pure water interface.
  ! If a pair of water molecules are h-bonded at time t, and both of them are located 
  ! in one of interfaces, then h(t)=1, else h(t)=0
  ! 0a) Date     Programmer     version
  !  ======      ==========     =======
  !  2020/9/18   Gang Huang     Original code
  !============
  ! Declaration
  !============
  USE parameter_shared,ONLY:n_grid
  USE atom_module
  USE count_time, ONLY:begin_time,end_time,rat
  USE types
  USE surf_traj
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  INTEGER,PARAMETER :: rk=8              
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
  CHARACTER(LEN=200) :: surf_filename
  CHARACTER(LEN=200) :: filename, pos_filename,list_oxygen_pairs
  CHARACTER(LEN=2) :: guest_atom, host_atom
  REAL(KIND=rk), DIMENSION(3) :: boxsize
  INTEGER :: criterion
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
  list_oxygen_pairs=""
  CALL system_clock(begin_time,rat) !To get the starting time
  
  !============================================
  ! To read the required controlling parameters
  !============================================
  CALL read_interface_input(boxsize,delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,&
           criterion,surf_filename) 
  !========================
  ! Sampling the trajectory
  !========================
  !CASE1: If we do not need to recenter, we can call sample_format2()
  !CALL sample_format2(pos_filename,nmo_start,nmo_end,nat,ns,n_samples)
  !CASE2: If we have to recenter, we call sample_and_recenter_format2()
  CALL load(pos_filename,nmo_start,nmo_end,nat,ns,n_samples,boxsize)
  ! After running the sample() subroutine, now we have a new sampled trajectory file 
  !(stored in atom_info), which is generally shorter than the original one.
  
  !====================
  !read surf trajectory
  !====================
  open(20,file=trim(surf_filename))
      CALL read_surf_traj(20,nmo_start,nmo_end,ns,n_grid,n_samples)
      !CALL read_surf_traj(20,0,999,1,n_grid,1000)
  close(20)

  CALL ghbond(filename,pos_filename,list_oxygen_pairs,nat) ! O-O pairs
  !Calculate n_HB(t) for pure water system. If the format of data is different, one may use another funtion, 
  !eg, ghbacf_n_pbc_format2(),which is very similar to this one.
  CALL ghbacf_interface_n_pbc_format2(boxsize,delta_t0,filename,pos_filename,list_oxygen_pairs, &
           n_samples,nat,ns,criterion)
  CALL system_clock(end_time,rat)
  WRITE(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 

END PROGRAM
