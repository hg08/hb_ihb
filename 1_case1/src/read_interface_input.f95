SUBROUTINE read_interface_input(boxsize,delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,&
               criterion,surf_filename,thickness) 
  !2018/12/27
  ! 0) Purpose:     
  ! 1a) The subroutine read_interface_input.f95 read the parameters from the input file. 
  ! One important input file is: the surface trajectory file. (And its name should be given to argument surf_filename.)
  ! The surface trajectory file is obtained from the program XXX 
  ! 0a) Date    Programmer     version
  ! 20181225    Gang Huang     3c
  ! Declaration
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  INTEGER,PARAMETER :: rk=8              

  !=========
  !Variables
  !=========
  real(KIND=rk), DIMENSION(3), INTENT(INOUT) :: boxsize  
  INTEGER, INTENT(INOUT) :: criterion
  REAL(KIND=rk), INTENT(INOUT) :: delta_t0  ! For reading data
  character(LEN=200), INTENT(INOUT) :: filename
  INTEGER :: i,iatom,ierr,imovie
  character(LEN=50) :: line 
  INTEGER, INTENT(INOUT) :: nat ! number of atoms
  INTEGER, INTENT(INOUT) :: nmo_start 
  INTEGER, INTENT(INOUT) :: nmo_end
  INTEGER, INTENT(INOUT) :: ns
  CHARACTER(LEN=200), INTENT(INOUT) :: pos_filename
  CHARACTER(LEN=200), INTENT(INOUT) :: surf_filename
  REAL(KIND=rk), INTENT(INOUT) :: thickness ! the thickness of the instantaneous interfaces

  !==============
  !Initialization
  !==============
  iatom=0; imovie=0; i=0

  !==================
  !read data in input
  !==================
  !WRITE(6,*)'What is the size of box (a,b,c):'
  READ(*,'(a)') line
  if (len_trim(line)==0) then
      !boxsize(1) = 5.316
      !boxsize(2) = 5.316
      !boxsize(3) = 5.316
      boxsize(1) = 5.32d0
      boxsize(2) = 5.32d0
      boxsize(3) = 5.32d0
  ELSE
    READ(line, *, iostat = ierr) boxsize(1), boxsize(2), boxsize(3)
    WRITE(*,*)"boxsize:", boxsize(1), boxsize(2), boxsize(3)
  ENDIF

  !WRITE(6,*)'What is the time step in the traj. file (ps): (Default: 0.d0005)'
  READ(*,'(a)') line
  IF (len_trim(line)==0) THEN
    delta_t0 = 0.0005d0
  ELSE
    READ(line, *, iostat = ierr) delta_t0
    WRITE(*,*) delta_t0
  ENDIF

  !WRITE(6,*)'What is the name of the system:'
  READ(*,'(a)') line
  IF (len_trim(line)==0) THEN
    filename = "8w_bulk_pbc"
  ELSE
    READ(line, *, iostat = ierr) filename
    WRITE(*,*) "filename: ",filename
  ENDIF

  !WRITE(6,*)'What is the name of the trajecotry file:'
  READ(*,'(a)') line
  IF (len_trim(line)==0) THEN
    pos_filename = "traj_pos_8w.xyz"
  ELSE
    READ(line, *, iostat = ierr) pos_filename
    WRITE(*,*) "pos_filename: ",pos_filename
  ENDIF

  !WRITE(6,*)'What is the initial step of the trajecotry:'
  READ(*,'(a)') line
  IF (len_trim(line)==0) THEN
    nmo_start = 0
  ELSE
    READ(line, *, iostat = ierr) nmo_start
    WRITE(*,*) "nmo_start: ", nmo_start
  ENDIF

  !WRITE(6,*)'What is the end step of the trajecotry:'
  READ(*,'(a)') line
  IF (len_trim(line)==0) THEN
    nmo_end = 0
  ELSE
    READ(line, *, iostat = ierr) nmo_end
    WRITE(*,*) "nmo_end: ", nmo_end
  ENDIF

  !WRITE(6,*)'What is the total number of atoms in the system:'
  READ(*,'(a)') line
  IF (len_trim(line)==0) THEN
    nat = 18
  ELSE
    READ(line, *, iostat = ierr) nat
    WRITE(*,*) "nat:", nat
  ENDIF

  !WRITE(6,*)'What is the time step for calculating CORRELATION:' 
  ! [ns*0.d0005] ps is the new time step for calculating correl func.
  READ(*,'(a)') line
  IF (len_trim(line)==0) THEN
    ns = 40
  ELSE
    READ(line, *, iostat = ierr) ns
    WRITE(*,*) "ns: ", ns
  ENDIF

  !WRITE(6,*)'What is the criterion of HB:' 
  READ(*,'(a)') line
  IF (len_trim(line)==0) THEN
    criterion = 1
  ELSE
    READ(line, *, iostat = ierr) criterion
    WRITE(*,*) "HB criterion: ", criterion
    WRITE(*,*) "criterion 1 denote ADH and 2 denotes AHD definition. "
  ENDIF  

  !WRITE(6,*)'What is the name of the surface trajectory file:' 
  READ(*,'(a)') line
  IF (len_trim(line)==0) THEN
    surf_filename = "surf_traj.dat"
  ELSE
    READ(line, *, iostat = ierr) surf_filename
    WRITE(*,*) "surf traj. name: ", surf_filename
  ENDIF  

  !WRITE(6,*)'What is the thickness of the interface you want to define: (Angstrom)' 
  READ(*,'(a)') line
  IF (len_trim(line)==0) THEN
    thickness = 3.d0
  ELSE
    READ(line, *, iostat = ierr) thickness
    !WRITE(*,*) "READ thickness: ", thickness
    !WRITE(*,*) "READ thickness: ", thickness
    !WRITE(*,*) "READ thickness: ", thickness
  ENDIF  

END SUBROUTINE read_interface_input
