SUBROUTINE read_interface_input_2(boxsize,delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,&
             surf_filename,thickness) 
  !2018/12/27
  ! 0) Purpose:     
  ! 1a) The subroutine read_interface_input_2.f95 read the parameters from the input file. 
  ! In read_interface_input_2.f95, criterion for HB definition is not needed.
  ! One important input file is: the surface trajectory file. (And its name should be given to argument surf_filename.)
  ! The surface trajectory file is obtained from the program XXX 
  ! 0a) Date    Programmer     version
  ! 20181225    Gang Huang     3c
  ! Declaration
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  integer,parameter :: rk=8              

  !=========
  !Variables
  !=========
  real(kind=rk), dimension(3), INTENT(INOUT) :: boxsize  
  !integer, INTENT(INOUT) :: criterion
  REAL(KIND=rk), INTENT(INOUT) :: delta_t0  ! For reading data
  character(LEN=200), INTENT(INOUT) :: filename
  INTEGER :: i,iatom,ierr,imovie
  character(LEN=200) :: line 
  INTEGER, INTENT(INOUT) :: nat ! number of atoms
  INTEGER, INTENT(INOUT) :: nmo_start 
  INTEGER, INTENT(INOUT) :: nmo_end
  INTEGER, INTENT(INOUT) :: ns
  character(LEN=200), INTENT(INOUT) :: pos_filename
  character(LEN=200), INTENT(INOUT) :: surf_filename
  REAL(KIND=rk), INTENT(INOUT) :: thickness ! the thickness of the instantaneous interfaces

  !==============
  !Initialization
  !==============
  iatom=0;imovie=0;i=0

  !==================
  !read data in input
  !==================
  write(6,*)'What is the size of box (a,b,c):'
  read(*,'(a)') line
  if (len_trim(line)==0) then
      boxsize(1) = 5.32d0
      boxsize(2) = 5.32d0
      boxsize(3) = 5.32d0
  else
    read(line, *, iostat = ierr) boxsize(1), boxsize(2), boxsize(3)
    write(*,*)"boxsize:", boxsize(1), boxsize(2), boxsize(3)
  endif

  write(6,*)'What is the time step in the traj. file (ps): (Default: 0.d0005)'
  read(*,'(a)') line
  if (len_trim(line)==0) then
    delta_t0 = 0.0005d0
  else
    read(line, *, iostat = ierr) delta_t0
    !write(*,*) delta_t0
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
    !write(*,*) "nat:", nat
  endif

  write(6,*)'What is the time step for calculating CORRELATION:' 
  ! [ns*0.d0005] ps is the new time step for calculating correl func.
  read(*,'(a)') line
  if (len_trim(line)==0) then
    ns = 40
  else
    read(line, *, iostat = ierr) ns
    !write(*,*) "ns: ", ns
  endif

  write(6,*)'What is the name of the surface trajectory file:' 
  read(*,'(a)') line
  if (len_trim(line)==0) then
    surf_filename = "surf_traj.dat"
  else
    read(line, *, iostat = ierr) surf_filename
    !write(*,*) "surf traj. name: ", surf_filename
  endif  

  write(6,*)'What is the thickness of the interface you want to define: (Angstrom)' 
  read(*,'(a)') line
  if (len_trim(line)==0) then
    thickness = 3.d0
  else
    read(line, *, iostat = ierr) thickness
    write(*,*) "READ thickness: ", thickness
  endif  

END SUBROUTINE read_interface_input_2 
