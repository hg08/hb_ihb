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
  if (len_trim(line)==0) then
      boxsize(1) = 5.316
      boxsize(2) = 5.316
      boxsize(3) = 5.316
  else
    read(line, *, iostat = ierr) boxsize(1), boxsize(2), boxsize(3)
    write(*,*) boxsize(1), boxsize(2), boxsize(3)
  endif
  write(6,*)'What is the time step in the traj. file (ps): (Default: 0.0005)'
  !read(5,*)delta_t0
  read(*,'(a)') line
  if (len_trim(line)==0) then
    delta_t0 = 0.0005
  else
    read(line, *, iostat = ierr) delta_t0
    write(*,*) delta_t0
  endif
  write(6,*)'What is the name of the system:'
  read(*,'(a)') line
  if (len_trim(line)==0) then
    filename = "8w_bulk_pbc"
  else
    read(line, *, iostat = ierr) filename
    write(*,*) "filename: ",filename
  endif
  write(6,*)'What is the name of the trajecotry file:'
  read(*,'(a)') line
  if (len_trim(line)==0) then
    pos_filename = "traj_pos_8w.xyz"
  else
    read(line, *, iostat = ierr) pos_filename
    write(*,*) "pos_filename: ",pos_filename
  endif
  write(6,*)'What is the initial step of the trajecotry:'
  read(*,'(a)') line
  if (len_trim(line)==0) then
    nmo_start = 0
  else
    read(line, *, iostat = ierr) nmo_start
    write(*,*) "nmo_start: ", nmo_start
  endif
  write(6,*)'What is the end step of the trajecotry:'
  read(*,'(a)') line
  if (len_trim(line)==0) then
    nmo_end = 0
  else
    read(line, *, iostat = ierr) nmo_end
    write(*,*) "nmo_end: ", nmo_end
  endif
  write(6,*)'What is the total number of atoms in the system:'
  read(*,'(a)') line
  if (len_trim(line)==0) then
    nat = 18
  else
    read(line, *, iostat = ierr) nat
    write(*,*) "nat:", nat
  endif
  write(6,*)'What is the time step for calculating CORRELATION:' 
  ! [ns*0.0005] ps is the new time step for calculating correl func.
  read(*,'(a)') line
  if (len_trim(line)==0) then
    ns = 10
  else
    read(line, *, iostat = ierr) ns
    write(*,*) "ns: ", ns
  endif
  write(6,*)'What is the criterion of HB:' 
  read(*,'(a)') line
  if (len_trim(line)==0) then
    criterion = 1
  else
    read(line, *, iostat = ierr) criterion
    write(*,*) "HB criterion: ", criterion
    write(*,*) "criterion 1 denote ADH and 2 denotes AHD definition. "
  endif  

  write(6,*)'What is the name of the surface trajectory file:' 
  read(*,'(a)') line
  if (len_trim(line)==0) then
    surf_filename = "surf_traj.dat"
  else
    read(line, *, iostat = ierr) surf_filename
    write(*,*) "surf traj. name: ", surf_filename
  endif  

  write(6,*)'What is the thickness of the interface you want to define: (Angstrom)' 
  read(*,'(a)') line
  if (len_trim(line)==0) then
    thickness = 3.0
  else
    read(line, *, iostat = ierr) thickness
    write(*,*) "thickness: ",  thickness
  endif  

END SUBROUTINE read_interface_input 
