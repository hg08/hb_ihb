SUBROUTINE read_interface_input(boxsize,delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,&
               criterion,surf_filename) 
  !2018/12/27
  ! 0) Purpose:     
  ! 0a) Date     Programmer     version
  !  20181225    Gang Huang     3c
  ! 1a) The subroutine read_interface_input.f95 read the parameters from the input file. 
  ! Declaration
  USE parameter_shared

  IMPLICIT NONE
  !==========
  !parameters
  !==========
  integer,parameter :: rk=8              
  INTEGER :: iatom,imovie,i
  character(LEN=50) :: line 
  real(kind=rk), dimension(3), INTENT(INOUT) :: boxsize  
  REAL(kind=rk), INTENT(INOUT) :: delta_t0  ! For reading data
  character(LEN=200), INTENT(INOUT) :: filename
  character(LEN=200), INTENT(INOUT) :: pos_filename
  character(LEN=200), INTENT(INOUT) :: surf_filename
  INTEGER, INTENT(INOUT) :: nmo_start 
  INTEGER, INTENT(INOUT) :: nmo_end
  INTEGER, INTENT(INOUT) :: nat ! number of atoms
  INTEGER, INTENT(INOUT) :: ns
  integer, INTENT(INOUT) :: criterion
  integer :: ierr

  !==============
  !Initialization
  !==============
  iatom = 0
  imovie =0
  i =0
  !==================
  !read data in input
  !==================
  write(6,*)'What is the size of box (a,b,c):'
  read(*,'(a)') line
  read(line, *, iostat = ierr) boxsize(1), boxsize(2), boxsize(3)
  write(*,*) boxsize(1), boxsize(2), boxsize(3)
  write(6,*)'What is the time step in the traj. file (ps): (Default: 0.0005)'

  read(*,'(a)') line
  read(line, *, iostat = ierr) delta_t0
  write(*,*) delta_t0
  write(6,*)'What is the name of the system:'

  read(*,'(a)') line
  read(line, *, iostat = ierr) filename
  write(*,*) "filename: ",filename
  write(6,*)'What is the name of the trajecotry file:'

  read(*,'(a)') line
  read(line, *, iostat = ierr) pos_filename
  write(*,*) "pos_filename: ",pos_filename
  write(6,*)'What is the initial step of the trajecotry:'

  read(*,'(a)') line
  read(line, *, iostat = ierr) nmo_start
  write(*,*) "nmo_start: ", nmo_start
  write(6,*)'What is the end step of the trajecotry:'
  
  read(*,'(a)') line
  read(line, *, iostat = ierr) nmo_end
  write(*,*) "nmo_end: ", nmo_end
  write(6,*)'What is the total number of atoms in the system:'

  read(*,'(a)') line
  read(line, *, iostat = ierr) nat
  write(*,*) "nat:", nat
  write(6,*)'What is the time step for calculating CORRELATION:' 

  read(*,'(a)') line
  read(line, *, iostat = ierr) ns
  write(*,*) "ns: ", ns
  write(6,*)'What is the criterion of HB:' 

  read(*,'(a)') line
  read(line, *, iostat = ierr) criterion
  write(*,*) "HB criterion: ", criterion
  write(*,*) "criterion 1 denote ADH and 2 denotes AHD definition. "

  write(6,*)'What is the name of the surface trajectory file:' 
  read(*,'(a)') line
  read(line, *, iostat = ierr) surf_filename
  write(*,*) "surf traj. name: ", surf_filename

  write(6,*)'What is the thickness of the interface you want to define: (Angstrom)' 
  read(*,'(a)') line
  read(line, *, iostat = ierr) thickness
  write(*,*) "thickness: ",  thickness

END SUBROUTINE read_interface_input 
