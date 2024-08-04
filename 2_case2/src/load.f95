SUBROUTINE load(pos_filename,nmo_start,nmo_end,nat,ns,n_samples,boxsize)
  !2018/12/27
  ! 0) Purpose:     
  ! 0a) The subroutine sample_and_recenter.f95 reduce the size of the trajectory and recenter each step of the trajectory. 
  ! 0b) Date     Programmer     version
  ! 2020-9-16   Gang Huang     1
  ! Declaration
  USE parameter_shared, ONLY: sampled_movie,sampled_time, &
      nb_divx, nb_divy, nb_divz, n_grid, divx, divy, divz, whish_size
  USE atom_module, ONLY: atom_info
  USE traj_format2
  IMPLICIT NONE
  !==========
  !Parameters
  !==========
  integer, parameter :: rk=8
  character(LEN=*), INTENT(IN) :: pos_filename
  INTEGER, INTENT(IN) :: nmo_start  
  INTEGER, INTENT(IN) :: nmo_end  
  INTEGER, INTENT(IN) :: nat ! number of atoms or nb_atoms
  INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
  INTEGER, INTENT(OUT) :: n_samples  !n_samples = INT(nmo/ns)
  real(kind=rk), dimension(3), INTENT(INOUT) :: boxsize
  !Local varables
  INTEGER :: iatom,imovie,i
  !INTEGER :: num_wat_pairs 
  INTEGER, PARAMETER :: n_axis = 2  ! n_axis=2 means z-axis is the normal axis.
  REAL(KIND=rk) :: sum_mass
  REAL(KIND=rk) :: center_pos 
  !Initialization
  iatom = 0
  imovie = 0
  i = 0 
  !num_wat_pairs = 0
  nb_divx = nint(boxsize(1)/whish_size) ! round the argument to the nearest integer.
  nb_divy = nint(boxsize(2)/whish_size) ! - 
  nb_divz = nint(boxsize(3)/whish_size) ! -
  divx = boxsize(1)/nb_divx
  divy = boxsize(2)/nb_divy
  divz = boxsize(3)/nb_divz
  n_grid = nb_divx * nb_divy

  n_samples = sampling_number(nmo_start,nmo_end,ns)
  allocate(sampled_movie(n_samples))
  allocate(sampled_time(n_samples))
  !=======================
  !read in trajectory file 
  !=======================
  open(10,file=trim(pos_filename))
  ! Now starting read data
  write(*,*) "Total number of atoms: ",nat
  CALL read_traj(10,nmo_start,nmo_end,ns,nat,n_samples) 
  close(10)
  write(6,*) 'End of trajectory reading.'

  deallocate(sampled_movie, sampled_time)
END SUBROUTINE load 
