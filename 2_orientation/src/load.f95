SUBROUTINE load(pos_filename,nmo_start,nmo_end,nat,ns,n_samples,boxsize, &
           sampled_pos_filename,sampled_movie,sampled_time, &
           nb_divx,nb_divy,nb_divz,n_grid,divx,divy,divz,whish_size,atom_info)
  !2018/12/27
  ! 0) Purpose:     
  ! 0a) The subroutine load.f95 loads the trajectory. 
  ! 0b) Date     Programmer     version
  ! 2020-9-16    Gang Huang     1
  ! 2024-8-5     Gang Huang     2
  ! Declaration
  USE module_ihb, ONLY: atom,read_traj
      
  IMPLICIT NONE
  !==========
  !Parameters
  !==========
  integer, parameter :: rk=8

  ! Varables
  character(LEN=200), INTENT(IN) :: pos_filename
  INTEGER, INTENT(IN) :: nmo_start  
  INTEGER, INTENT(IN) :: nmo_end  
  INTEGER, INTENT(IN) :: nat ! number of atoms or nb_atoms
  INTEGER, INTENT(IN) :: ns ! Get one sample from the trajectory every ns step.
  INTEGER, INTENT(IN) :: n_samples !n_samples = INT(nmo/ns)
  real(kind=rk), dimension(3), INTENT(INOUT) :: boxsize
  REAL(KIND=rk), INTENT(IN) :: whish_size ! Angstrom
  character(LEN=200), INTENT(INOUT):: sampled_pos_filename
  INTEGER, DIMENSION(n_samples),INTENT(INOUT) :: sampled_movie
  REAL(KIND=rk), DIMENSION(n_samples), INTENT(INOUT) :: sampled_time
  INTEGER, INTENT(INOUT) :: nb_divx, nb_divy, nb_divz, n_grid 
  REAL(KIND=rk), INTENT(INOUT) :: divx, divy, divz
  TYPE(atom), DIMENSION(nat,n_samples),INTENT(INOUT) :: atom_info
  ! Local varables
  INTEGER :: iatom,imovie,i
  INTEGER :: num_wat_pairs 
  INTEGER, PARAMETER :: n_axis = 2  ! n_axis=2 means z-axis is the normal axis.
  REAL(KIND=rk) :: sum_mass
  REAL(KIND=rk) :: center_pos 
  
  ! Initialization
  iatom = 0; imovie = 0; i = 0; num_wat_pairs = 0
  nb_divx = nint(boxsize(1)/whish_size) ! round the argument to the nearest integer.
  nb_divy = nint(boxsize(2)/whish_size) ! round the argument to the nearest integer.
  nb_divz = nint(boxsize(3)/whish_size) ! round the argument to the nearest integer.
  divx = boxsize(1)/REAL(nb_divx,rk)
  divy = boxsize(2)/REAL(nb_divy,rk)
  divz = boxsize(3)/REAL(nb_divz,rk)
  n_grid = nb_divx * nb_divy
  num_wat_pairs = (nat/3)*(nat/3-1)/2

  !=======================
  !read in trajectory file 
  !=======================
  open(10,file=trim(pos_filename))
  ! Now starting read data
  CALL read_traj(10,nmo_start,nmo_end,ns,nat,n_samples,sampled_movie,sampled_time,atom_info) 
  close(10)
  write(6,*) 'End of trajectory reading.'

  !=============
  !write in file
  !=============
  sampled_pos_filename = 'recentered_traj_pos_sampled.xyz'
  open(10,file=sampled_pos_filename)

  step: DO i=1,n_samples
    sum_mass = 0.d0
    center_pos = 0.d0

    write (10,'(I8)') nat
    WRITE(10,100) ' i = ',i-1,', time = ',sampled_time(i)
    100 FORMAT (A5,I8,A9,F12.3)

    ! Format for writing in recentered 
    131 FORMAT (A4,3F20.10)
     
    DO iatom = 1, nat
        WRITE(10,131) TRIM(atom_info(iatom, i)%atom_name), &
        atom_info(iatom,i)%coord(1), &
        atom_info(iatom,i)%coord(2), &
        atom_info(iatom,i)%coord(3)
    ENDDO
  ENDDO step

  write(6,*)'Sampled trajectory is written: ', sampled_pos_filename
  close(10)

END SUBROUTINE 
