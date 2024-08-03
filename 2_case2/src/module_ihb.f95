MODULE atom_module
! Purpose:
!   To define the derived data type for atom
!
IMPLICIT NONE
TYPE :: atom
  CHARACTER(LEN=2) :: atom_name
  INTEGER :: atom_id 
  INTEGER :: host_id  ! For O atom in water, host_id = atom_id
  REAL :: mass
  REAL, DIMENSION(3) :: coord 
END TYPE atom

! The array atom_info can be shared by subroutines  
TYPE(atom), ALLOCATABLE, DIMENSION(:,:) :: atom_info
END MODULE atom_module
MODULE count_time
! Purpose:
!   To declear data for calculating the total time cosumme in the process.
!
IMPLICIT NONE
integer :: begin_time,end_time,rat
END MODULE count_time
MODULE parameter_shared
!
! Purpose:
!   To declare data to share between routines.
IMPLICIT NONE
SAVE 
character(LEN=200) :: sampled_pos_filename
INTEGER, ALLOCATABLE, DIMENSION(:) :: sampled_movie
REAL, ALLOCATABLE, DIMENSION(:) :: sampled_time, sampled_energy

! For griding
REAL, PARAMETER :: whish_size=0.5 ! Angstrom
INTEGER :: nb_divx, nb_divy, nb_divz, n_grid 
REAL :: divx, divy, divz

! paramters of the interface
REAL :: thickness ! the thickness of the instantaneous interfaces

END MODULE parameter_shared
MODULE surf_module
! Purpose:
!   To define the derived data type for isosurfaces of an interfacial system
IMPLICIT NONE
TYPE :: surf
  INTEGER :: grid_index_x, grid_index_y 
  REAL(KIND=8), DIMENSION(2) :: coord ! upper and lower location of the isosurface two real values.
END TYPE surf

! The array surf_info can be shared by subroutines  
TYPE(surf), ALLOCATABLE, DIMENSION(:,:) :: surf_info
END MODULE surf_module
MODULE surf_traj
!
! Purpose: 
! To declear data related to the simulation and traj.
!
!========================================
! data         author         version
! ========     =======        ==========
! 2020-9-18    Gang Huang     Original
!========================================
!
IMPLICIT NONE

CONTAINS

  SUBROUTINE read_surf_traj(indx,nmo_start,nmo_end,ns,n_grid,n_samples)
    ! Purpose:
    ! To read info from the isosurface trajectory file (format: surf_traj_*.xyz)
    USE surf_module, ONLY: surf_info
    INTEGER :: ierror, i_sample, i_grid, ix, iy
    INTEGER, INTENT(IN) :: n_grid  ! KEEP
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns) KEEP
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER, INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    CHARACTER(LEN=4) :: head_char,e
    INTEGER :: y

    ALLOCATE(surf_info(n_grid,n_samples))
    i_sample = 0
    IF (nmo_start > 0) THEN
          DO i_sample = 0, nmo_start ! Skip 
               read (indx,*) 
               do i_grid = 1, n_grid
                  read (indx,*) 
               enddo
          END DO
    ENDIF
    DO i_sample = 0, n_samples - 1
          read(indx, *) head_char, e, y 
          !write(*,*), head_char, e, y
          do i_grid = 1, n_grid
              read(indx, *, IOSTAT=ierror) ix, iy, surf_info(1,i_grid,i_sample+1), surf_info(2,i_grid,i_sample+1)
              !write(*,*) 'ix iy Surf1 Surf2', ix, iy, surf_info(1,i_grid,i_sample+1), surf_info(2,i_grid,i_sample+1)
          enddo 
    END DO 
    close(indx)
  END SUBROUTINE read_surf_traj
  !SUBROUTINE read_surf_traj(indx,nmo_start,nmo_end,ns,n_grid,n_samples)
  !  ! Purpose:
  !  ! To read info from the isosurface trajectory file (format: surf_traj_*.xyz)
  !  USE surf_module, ONLY: surf_info
  !  INTEGER :: i_sample, i_grid
  !  INTEGER, INTENT(IN) :: n_grid  ! KEEP
  !  INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns) KEEP
  !  INTEGER, INTENT(IN) :: indx
  !  INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
  !  INTEGER, INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
  !  CHARACTER(LEN=4) :: head_char
  !  INTEGER :: y

  !  allocate(surf_info(n_grid,n_samples))
  !  i_sample = 1
  !  write(*,*) "read_surf_traj(): New total time steps (n_samples):", n_samples
  !  DO WHILE (i_sample < n_samples+1) ! +1 means i_sample can take the value of n_samples 
  !      read(indx, '(A4)') head_char
  !      PRE_CHECK:IF (head_char=="i = ") THEN
  !          BACKSPACE(UNIT=indx) ! Because I am not able to read other lines with the format '(A4,I8)', and have not find any good way, so I try to read it in '(A4)' first 
  !          read(indx, '(A4,I5)') head_char, y
  !          !read(indx, '(1X,A4,I5)') head_char, y
  !          CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-(nmo_start-1),ns) == 0) THEN
  !              !-------------------------------------------------------------------------------------------------------
  !              !NOTE: if use ' i = ', instead of 'i = ', it will be wrong!
  !              !IF (head_char==' i = ' .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-(nmo_start-1),ns) == 1) THEN
  !              !-------------------------------------------------------------------------------------------------------
  !              WRITE(*,*)"read_traj():", head_char, y
  !              BACKSPACE(UNIT=indx) ! Because we have to read the whole line with ' i = ' line.
  !              READ(indx,*) ! skip one line in the unit=indx
  !              131 FORMAT (11X,2F13.6)
  !              inner: do i_grid= 1, n_grid
  !                read (indx,131) surf_info(i_grid,i_sample)%coord(1), & 
  !                  surf_info(i_grid,i_sample)%coord(2)
  !                  WRITE (*,131) surf_info(i_grid,i_sample)%coord(1), &
  !                  surf_info(i_grid, i_sample)%coord(2)
  !              enddo inner
  !              i_sample = i_sample + 1 ! The position is important. It must be located before ENDIF MATCH
  !          ENDIF CHECK_HEAD
  !      ENDIF PRE_CHECK
  !  END DO
  !END SUBROUTINE read_surf_traj
END MODULE surf_traj

      MODULE tools 
      !2020/2/15
      !===================================================
      ! The functions 
      !===================================================      
      implicit none

      integer,private,parameter :: rk=4  

      CONTAINS
 
      REAL FUNCTION direct_dist2(u1,v1,w1,u2,v2,w2)
          real(KIND=rk), INTENT(IN) :: u1,v1,w1,u2,v2,w2
          direct_dist2 = (u2-u1)**2 + (v2-v1)**2 + (w2-w1)**2
      END FUNCTION direct_dist2
      
      REAL FUNCTION distance2(r1,r2,boxsize)
          real(KIND=rk), DIMENSION(3), INTENT(IN) :: r1,r2
          real(KIND=rk), DIMENSION(3), INTENT(IN) :: boxsize
          REAL(KIND=rk) :: dx,dy,dz
          dx = r1(1) - r2(1)
          if (abs(dx) > boxsize(1)*0.5) then
              dx = boxsize(1) - dx
          endif
          dy = r1(2) - r2(2)
          if (abs(dy) > boxsize(2)*0.5) then
              dy = boxsize(2) - dy
          endif
          dz = r1(3) - r2(3)
          if (abs(dz) > boxsize(3)*0.5) then
              dz = boxsize(3) - dz
          endif
          distance2 = dx**2 + dy**2 + dz**2
      END FUNCTION distance2
     
      !REAL FUNCTION dist2(u1,v1,w1,u2,v2,w2,a,b,c)
      !    real(KIND=rk),INTENT(INOUT) :: u1,v1,w1,u2,v2,w2,a,b,c
      !    logical :: A1,A2,A3,B1,B2,B3,C1,C2,C3
      !    A1 = (abs(u1-u2) > a/2 .AND. u1 > u2)
      !    A2 = (abs(u1-u2) > a/2 .AND. u2 > u1)
      !    A3 = (abs(u1-u2) < a/2)
      !    B1 = (abs(v1-v2) > b/2 .AND. v1 > v2)
      !    B2 = (abs(v1-v2) > b/2 .AND. v2 > v1)
      !    B3 = (abs(v1-v2) < b/2)
      !    C1 = (abs(w1-w2) > c/2 .AND. w1 > w2)
      !    C2 = (abs(w1-w2) > c/2 .AND. w2 > w1)
      !    C3 = (abs(w1-w2) < c/2)
      !    if (A3 .and. B3 .and. C3) then
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A2 .and. B3 .and. C3) then
      !        u2 = u2-a
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A1 .and. B3 .and. C3) then
      !        u1 = u1-a
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A3 .and. B1 .and. C3) then
      !        v1 = v1 - b
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A3 .and. B3 .and. C2) then
      !        w2 = w2 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A3 .and. B2 .and. C3) then
      !        v2 = v2 - b
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A3 .and. B3 .and. C1) then
      !        w1 = w1 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A1 .and. B1 .and. C3) then
      !        u1 = u1 - a
      !        v1 = v1 - b
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A1 .and. B2 .and. C3) then
      !        u1 = u1 - a
      !        v2 = v2 - b
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A1 .and. B3 .and. C1) then
      !        u1 = u1 - a
      !        w1 = w1 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A1 .and. B3 .and. C2) then
      !        u1 = u1 - a
      !        w2 = w2 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A2 .and. B1 .and. C3) then
      !        u2 = u2 - a
      !        v1 = v1 - b
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A2 .and. B2 .and. C3) then
      !        u2 = u2 - a
      !        v2 = v2 - b
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A2 .and. B3 .and. C1) then
      !        u2 = u2 - a
      !        w1 = w1 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A2 .and. B3 .and. C2) then
      !        u2 = u2 - a
      !        w2 = w2 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A3 .and. B1 .and. C1) then
      !        v1 = v1 - b
      !        w1 = w1 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A3 .and. B1 .and. C2) then
      !        v1 = v1 - b
      !        w2 = w2 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A3 .and. B2 .and. C1) then
      !        v2 = v2 - b
      !        w1 = w1 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A3 .and. B2 .and. C2) then
      !        v2 = v2 - b
      !        w2 = w2 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A1 .and. B1 .and. C1) then
      !        u1 = u1 - a
      !        v1 = v1 - b
      !        w1 = w1 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A1 .and. B1 .and. C2) then
      !        u1 = u1 - a
      !        v1 = v1 - b
      !        w2 = w2 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A1 .and. B2 .and. C1) then
      !        u1 = u1 - a
      !        v2 = v2 - b
      !        w1 = w1 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A1 .and. B2 .and. C2) then
      !        u1 = u1 - a
      !        v2 = v2 - b
      !        w2 = w2 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A2 .and. B1 .and. C1) then
      !        u2 = u2 - a
      !        v1 = v1 - b
      !        w1 = w1 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A2 .and. B1 .and. C2) then
      !        u2 = u2 - a
      !        v1 = v1 - b
      !        w2 = w2 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A2 .and. B2 .and. C1) then
      !        u2 = u2 - a
      !        v2 = v2 - b
      !        w1 = w1 - c
      !        dist2= direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    elseif (A2 .and. B2 .and. C2) then
      !        u2 = u2 - a
      !        v2 = v2 - b
      !        w2 = w2 - c
      !        dist2 = direct_dist2(u1,v1,w1,u2,v2,w2) 
      !    endif
      !END FUNCTION dist2

      REAL FUNCTION diff_axis(u1,u2,h)
          ! u2 is used as origin
          real(KIND=rk),INTENT(IN) :: u1, u2
          real(KIND=rk),INTENT(IN) :: h
          real(KIND=rk) :: du 
          du = u1-u2
          if (abs(du) > 0.5*h) then
              du = h - du 
          endif
          diff_axis = du
      END FUNCTION diff_axis
      
      REAL FUNCTION pmADH(r1,r2,r3,boxsize)
          real(KIND=rk),DIMENSION(3), INTENT(IN) :: r1,r2,r3
          real(KIND=rk),DIMENSION(3), INTENT(IN) :: boxsize 
          pmADH= diff_axis(r3(1),r2(1),boxsize(1))*           &
                 diff_axis(r1(1),r2(1),boxsize(1))+      &
                 diff_axis(r3(2),r2(2),boxsize(2))*      &
                 diff_axis(r1(2),r2(2),boxsize(2))+      &
                 diff_axis(r3(3),r2(3),boxsize(3))*      &
                 diff_axis(r1(3),r2(3),boxsize(3))       ! pm: point multiplication. 
      END FUNCTION pmADH

      REAL FUNCTION pmAHD(r1,r2,r3,boxsize)
          real(KIND=rk),DIMENSION(3), INTENT(IN) :: r1,r2,r3
          real(KIND=rk),DIMENSION(3), INTENT(IN) :: boxsize 
          pmAHD= diff_axis(r1(1),r3(1),boxsize(1))*           &
                 diff_axis(r2(1),r3(1),boxsize(1))+      &
                 diff_axis(r1(2),r3(2),boxsize(2))*      &
                 diff_axis(r2(2),r3(2),boxsize(2))+      &
                 diff_axis(r1(3),r3(3),boxsize(3))*      &
                 diff_axis(r2(3),r3(3),boxsize(3))       ! pm: point multiplication. 
      END FUNCTION pmAHD 
      
      REAL FUNCTION diff_axis_v1(u1,u2,a)
          logical :: A1,A2,A3
          real(KIND=rk),INTENT(INOUT) :: u1, u2
          real(KIND=rk),INTENT(IN) :: a
          A1 = (abs(u1-u2) > a/2 .AND. u1 > u2)
          A2 = (abs(u1-u2) > a/2 .AND. u2 > u1)
          A3 = (abs(u1-u2) < a/2)
          if (A3) then
              diff_axis_v1= u1 - u2
          elseif (A1) then
              u1 = u1 - a
              diff_axis_v1= u1-u2
          elseif (A2) then
              u2 = u2 - a
              diff_axis_v1= u1 - u2
          endif
      END FUNCTION diff_axis_v1
     
      INTEGER FUNCTION get_total_numer_of_lines(file_)
          INTEGER :: nlines, io
          CHARACTER (len=*),INTENT(IN) :: file_

          !Initialization
          nlines = 0 
          
          OPEN (12, file = file_, status="old")
          DO
              READ(12,*,iostat=io)
              IF (io/=0) EXIT
              nlines = nlines + 1
          END DO
          CLOSE (12)

          get_total_numer_of_lines = nlines
      END FUNCTION 

      FUNCTION hydrogen_ndx_list(ndx_oxygen_1, & 
                 ndx_oxygen_2, pos_filename, natoms,boxsize) 
      implicit none
      !============
      ! 2020/07
      !========================
      !Parameters and variables
      !========================
      INTEGER, INTENT(IN) :: ndx_oxygen_1, ndx_oxygen_2
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      real(KIND=rk), dimension(3), INTENT(IN) :: boxsize
      
      !Local 
      real(KIND=rk), dimension(3) :: r_a, r_b
      integer,parameter :: rk=4              
      real,parameter :: r_ohc=1.21   ! rOH (1.1**2)
      real           :: r ! distance between O and H.
      integer    :: i,j,nmovie,iatom,& 
                    m1,m2,m3,i_H,&
                    i1,i2,ii,jj,i_OW, num
      INTEGER, DIMENSION(4) :: hydrogen_ndx_list
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:)          :: ndx_OW,ndx_H
      ! Initialization
      i1=0; i2=0
      
      !==================
      !Read data in input
      !==================
      nmovie=1!Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      ! To obtain the # of H and O atoms
      i_H = get_number_of_hydrogen(pos_filename, natoms) ! By calling functions, we decouple different functions
      i_OW = get_number_of_oxygen(pos_filename, natoms)! By calling functions, we decouple different functions 
      !=======================================================
      ! Calculate the indices of O (H) atoms in water molecules
      !=======================================================
      allocate(ndx_OW(i_OW))    ! this should be put after i_OW is defined
      allocate(ndx_H(i_H))      ! this should be put after i_H is defined
      i=0; ii=0; jj=0
      OPEN(10,file=trim(pos_filename))
          REWIND(10)     
          read(10,*) !Neglect data of this line. Be careful, the two 'read' lines must be located before the 'do iatom=1,natoms' loop, otherwise, there will be an error.
          read(10,*)                  
          do iatom=1,natoms
              read (10,*)atom_type(iatom),x(iatom,nmovie),& 
                  y(iatom,nmovie),z(iatom,nmovie)
              if (trim(atom_type(iatom)) .eq. 'O')then
                     i=i+1      
                     ndx_OW(i)=iatom
              elseif(trim(atom_type(iatom)) .eq. 'H') then
                     ii=ii+1
                     ndx_H(ii)=iatom
              else
              endif
          enddo 
      CLOSE(10)
      deallocate(atom_type)
      !====================
      !Producing the H list      
      !====================
      m1=ndx_oxygen_1
      m2=ndx_oxygen_2
      num = 0 ! For counting
      ! Case one
      do j =1, 1 ! Consider one step
          do ii=1,i_H 
              m3=ndx_H(ii)
              ! I use pbc to consider the distance r
              r_a= (/x(m2,j),y(m2,j),z(m2,j)/)  ! Coordinate of O2
              r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
              r=distance2(r_a,r_b,boxsize)  ! squared distance between O2 and H
              if (r<r_ohc) then  ! If the H is bound to O2
                  num = num + 1
                  hydrogen_ndx_list(num) = m3  ! Add the index of the H atom into the H list for the water pair (O1, O2)
              endif
          enddo
      enddo
      ! Case two
      do j =1, 1
          do ii=1,i_H
              m3=ndx_H(ii)
              ! I use pbc to consider the distance r
              r_a= (/x(m1,j),y(m1,j),z(m1,j)/) ! Coordinate of O1
              r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
              r=distance2(r_a,r_b,boxsize) ! squared distance between O2 and H
              if (r<r_ohc) then ! If the H is bound to O2
                  num = num + 1
                  hydrogen_ndx_list(num) = m3 ! Add the index of the H atom into the H list for the water pair (O1, O2)
              endif
          enddo
      enddo
      deallocate(ndx_OW,ndx_H,x,y,z)
      !=============================
      END FUNCTION hydrogen_ndx_list

      FUNCTION hydrogen_ndx_list_XO(ndx_X, & 
               ndx_oxygen_2,pos_filename, natoms,boxsize) 
      implicit none
      !============
      ! In this function, we use ndx_X replace ndx_oxygen_1. 
      ! ndx_X is not bonded with any Hydrogen atoms
      !========================
      !Parameters and variables
      !========================
      INTEGER, INTENT(IN) :: ndx_X, ndx_oxygen_2
      character(LEN=200),INTENT(IN) :: pos_filename      ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      real(KIND=rk), dimension(3), INTENT(IN) :: boxsize
      
      !Local 
      real(KIND=rk), dimension(3) :: r_a, r_b
      integer,parameter :: rk=4              
      real,parameter :: r_ohc=1.21   ! rOH (1.1**2)
      real :: r ! distance between O and H.
      integer :: i,j,nmovie,iatom,& 
              m1,m2,m3,i_H,&
              i1,i2,ii,jj,i_O,num
      INTEGER, DIMENSION(2) :: hydrogen_ndx_list_XO ! the length is 2, because there are 2 Hydrogen are bonded to ndx_oxygen_2 but no Hydrogen is bonded to ndx_X
      real,allocatable,dimension (:,:) :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:) :: ndx_O,ndx_H
      ! Initialization
      i1=0; i2=0
      m1=0; m2=0; m3=0
      
      !==================
      !Read data in input
      !==================
      nmovie=1!Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      ! To obtain the # of H and O atoms
      i_H = get_number_of_hydrogen(pos_filename, natoms) ! By calling functions, we decouple different functions
      i_O = get_number_of_oxygen(pos_filename, natoms)
      !write(*,*) "natoms: ", natoms
      !write(*,*) "i_O: ", i_O
      !write(*,*) "i_H: ", i_H
      !=======================================================
      ! Calculate the indices of O (H) atoms in water molecules
      !=======================================================
      allocate(ndx_O(i_O)) ! this should be put after i_O is defined
      allocate(ndx_H(i_H)) ! this should be put after i_H is defined
      !allocate(ndx_I(i_I)) ! this should be put after i_I is defined
      i=0; ii=0; jj=0
      write(*,*) "pos_filenae:",pos_filename
      open(10,file=trim(pos_filename))
      REWIND(10)     
      read(10,*) !Neglect data of this line. Be careful, the two 'read' lines must be located before the 'do iatom=1,natoms' loop, otherwise, there will be an error.
      read(10,*)                  
      do iatom=1,natoms
          read (10,*)atom_type(iatom),x(iatom,nmovie),& 
                   y(iatom,nmovie),z(iatom,nmovie)
          if (trim(atom_type(iatom)) .eq. 'O') then
                 ndx_O(i+1)=iatom
                 i=i+1      
                 write(*,*) "iatom:(test O) ", iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ndx_H(ii+1)=iatom
                 ii=ii+1
                 write(*,*) "iatom:(test) ", iatom
          else
          endif
      enddo 
      close(10)
      deallocate(atom_type)
      !====================
      !Producing the H list      
      !====================
      m1=ndx_X
      m2=ndx_oxygen_2
      num = 0 ! For counting
      ! Case one
      do j =1, 1 ! Consider one step
          do ii=1,i_H 
              m3=ndx_H(ii)
              write(*,*) "i_H: ", i_H
              write(*,*) "m3: ", m3
              ! I use pbc to consider the distance r
              r_a= (/x(m2,j),y(m2,j),z(m2,j)/)  ! Coordinate of O2
              r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
              r=distance2(r_a,r_b,boxsize)  ! squared distance between O2 and H
              if (r<r_ohc) then  ! If the H is bound to O2
                  num = num + 1
                  hydrogen_ndx_list_XO(num) = m3  ! Add the index of the H atom into the H list for the water pair (O1, O2)
              endif
          enddo
      enddo
      ! Case two is canceled
      deallocate(ndx_O,ndx_H,x,y,z)
      !=============================
      END FUNCTION hydrogen_ndx_list_XO

     ! FUNCTION hydrogen_ndx_list_between_nitrate_water(ndx_oxygen_1, & 
     !            ndx_oxygen_2,pos_filename, natoms,boxsize) 
     ! implicit none
     ! !========================
     ! !Parameters and variables
     ! !========================
     ! INTEGER, INTENT(IN) :: ndx_oxygen_1, ndx_oxygen_2
     ! character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
     ! INTEGER, INTENT(IN) :: natoms
     ! real(KIND=rk), dimension(3), INTENT(IN) :: boxsize
     ! 
     ! !Local 
     ! real(KIND=rk), dimension(3) :: r_a, r_b
     ! integer,parameter :: rk=4              
     ! real,parameter :: r_ohc=1.21   ! rOH (1.1**2)
     ! real           :: r ! distance between O and H.
     ! integer    :: i,j,nmovie,iatom,& 
     !               imovie,m1,m2,m3,i_H,&
     !               i1,i2,ii,jj,i_O,i_OW,num
     ! INTEGER, DIMENSION(4) :: hydrogen_ndx_list
     ! real,allocatable,dimension (:,:)           :: x,y,z
     ! character(LEN=3),allocatable,dimension (:) :: atom_type
     ! integer,allocatable,dimension (:)          :: ndx_OW,ndx_H
     ! ! Initialization
     ! i2 = 0; i1=0
     ! 
     ! !==================
     ! !Read data in input
     ! !==================
     ! nmovie=1!Number of movie steps reqired.
     ! allocate(atom_type(natoms))
     ! allocate(x(natoms,nmovie))
     ! allocate(y(natoms,nmovie))
     ! allocate(z(natoms,nmovie))
     ! !=======================
     ! !read in trajectory file 
     ! !=======================
     ! open(12,file=trim(pos_filename))
     ! REWIND(12)     
     ! do imovie=1,nmovie  !Note that nmovie=1
     !     read(12,*)!Neglect data of this line
     !     read(12,*)                  
     !     i=0
     !     ii=0
     !     do iatom= 1,natoms
     !         read (12,*)atom_type(iatom),x(iatom,imovie),& 
     !                  y(iatom,imovie),z(iatom,imovie)
     !         if (trim(atom_type(iatom)) .eq. 'O') then
     !             i=i+1
     !         elseif(trim(atom_type(iatom)) .eq. 'H') then
     !             ii=ii+1
     !         endif
     !     enddo
     ! enddo
     ! !i_O=i
     ! i_H=ii
     ! i_O = get_number_of_oxygens(pos_filename, natoms)
     ! !write(6,*)'Number of OW','    Number of HW'
     ! !write(6,*)i_OW,i_H
     ! close(12)
     ! !=======================================================
     ! ! Calculate the number of O (H) atoms in water molecules
     ! !=======================================================
     ! allocate(ndx_OW(i_O))    ! this should be put after i_O is defined 
     ! allocate(ndx_H(i_H))     ! this should be put after i_H is defined
     ! i=0; ii=0; jj=0
     ! do iatom=1,natoms
     !     if (trim(atom_type(iatom)) .eq. 'O')then
     !            i=i+1      
     !            ndx_OW(i)=iatom
     !     elseif(trim(atom_type(iatom)) .eq. 'H') then
     !            ii=ii+1
     !            ndx_H(ii)=iatom
     !     else
     !     endif
     ! enddo 
     ! deallocate(atom_type)
     ! !====================
     ! !Producing the H list      
     ! !====================
     ! m1=ndx_oxygen_1
     ! m2=ndx_oxygen_2
     ! num = 0 ! For counting
     ! ! Case one
     ! do j =1, 1 ! Consider one step
     !     do ii=1,i_H 
     !         m3=ndx_H(ii)
     !         ! I use pbc to consider the distance r
     !         r_a= (/x(m2,j),y(m2,j),z(m2,j)/) ! Coordinate of O2
     !         r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
     !         r=distance2(r_a,r_b,boxsize) ! squared distance between O2 and H
     !         if (r<r_ohc) then  ! If the H is bound to O2
     !             num = num + 1
     !             hydrogen_ndx_list(num) = m3  ! Add the index of the H atom into the H list for the water pair (O1, O2)
     !         endif
     !     enddo
     ! enddo
     ! ! Case two
     ! do j =1, 1
     !     do ii=1,i_H
     !         m3=ndx_H(ii)
     !         ! I use pbc to consider the distance r
     !         r_a= (/x(m1,j),y(m1,j),z(m1,j)/) ! Coordinate of O1
     !         r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
     !         r=distance2(r_a,r_b,boxsize) ! squared distance between O2 and H
     !         if (r<r_ohc) then ! If the H is bound to O2
     !             num = num + 1
     !             hydrogen_ndx_list(num) = m3 ! Add the index of the H atom into the H list for the water pair (O1, O2)
     !         endif
     !     enddo
     ! enddo
     ! deallocate(ndx_OW,ndx_H,x,y,z)
     ! !=============================
     ! END FUNCTION hydrogen_ndx_list_between_nitrate_water
      !TODO: define get_number_of_oxygens_in_nitrate()
      !TODO: define get_number_of_oxygens_in_water()
      !TODO: define get_number_of_iodine() Done
      FUNCTION get_number_of_oxygen(pos_filename, natoms) 
      implicit none
      !Parameters and variables
      !========================
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER :: get_number_of_oxygen 
      !Local 
      integer    :: i,nmovie,iatom,imovie
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1!Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          read(10,*)!Neglect data of this line
          read(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O') then
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_oxygen = i
      close(10)
      deallocate(atom_type,x,y,z)
      !=============================
      END FUNCTION get_number_of_oxygen

      FUNCTION get_number_of_hydrogen(pos_filename, natoms) 
      implicit none
      !Parameters and variables
      !========================
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER :: get_number_of_hydrogen
      !Local 
      integer    :: i,nmovie,iatom,imovie
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1!Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          read(10,*)!Neglect data of this line
          read(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'H') then
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_hydrogen = i
      close(10)
      deallocate(atom_type,x,y,z)
      !=============================
      END FUNCTION get_number_of_hydrogen

      FUNCTION get_number_of_iodine(pos_filename, natoms) 
      implicit none
      !Parameters and variables
      !========================
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER :: get_number_of_iodine 
      !Local 
      integer    :: i,nmovie,iatom,imovie
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1!Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          read(10,*)!Neglect data of this line
          read(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'I') then
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_iodine = i
      close(10)
      deallocate(atom_type,x,y,z)
      !=============================
      END FUNCTION get_number_of_iodine

      INTEGER FUNCTION grid_index(x,y,divx,divy,n_divy)
          ! transfer the coordinates (x,y) to grid_index, which is an integer
          REAL, INTENT(IN) :: x,y
          REAL, INTENT(IN) :: divx, divy
          INTEGER, INTENT(IN) :: n_divy
          INTEGER, DIMENSION(2) :: ind
          
          !Initialization
          ind = 0.0

          ind(1) = FLOOR(x/divx) 
          ind(2) = FLOOR(y/divy)
          
          grid_index = ind(1) * n_divy + ind(2)+1
    
      END FUNCTION grid_index

      LOGICAL FUNCTION pair_in_surf1(surf1_mol1,z1,surf1_mol2,z2,thickness)
          REAL(KIND=8), INTENT(IN) :: surf1_mol1,surf1_mol2
          REAL, INTENT(IN) :: z1,z2,thickness
          LOGICAL :: mol1_in_surf1, mol2_in_surf1

          mol1_in_surf1 = surf1_mol1 + & 
              thickness > z1
          mol2_in_surf1 = surf1_mol2 + &
              thickness > z2
          pair_in_surf1 = mol1_in_surf1 .and. mol2_in_surf1

      END FUNCTION pair_in_surf1    

      LOGICAL FUNCTION pair_in_surf2(surf2_mol1,z1,surf2_mol2,z2,thickness)
          REAL(KIND=8), INTENT(IN) :: surf2_mol1,surf2_mol2
          REAL, INTENT(IN) :: z1,z2,thickness
          LOGICAL :: mol1_in_surf2, mol2_in_surf2

          mol1_in_surf2 = surf2_mol1 - & 
              thickness < z1
          mol2_in_surf2 = surf2_mol2 - &
              thickness < z2
          pair_in_surf2 = mol1_in_surf2 .and. mol2_in_surf2

      END FUNCTION pair_in_surf2    

      LOGICAL FUNCTION mol_in_surf1(surf1_mol,z1,thickness)
          !Check if a molecule is in the surf1 ( the lower layer)
          REAL(KIND=8), INTENT(IN) :: surf1_mol
          REAL, INTENT(IN) :: z1,thickness
          !LOGICAL :: mol_in_surf1

          mol_in_surf1 = surf1_mol + & 
              thickness > z1

      END FUNCTION mol_in_surf1    

      LOGICAL FUNCTION mol_in_surf2(surf2_mol,z2,thickness)
          !Check if a molecule is in the surf2 ( the upper layer)
          REAL(KIND=8), INTENT(IN) :: surf2_mol
          REAL, INTENT(IN) :: z2,thickness
          !LOGICAL :: mol_in_surf2

          mol_in_surf2 = surf2_mol - & 
              thickness < z2

      END FUNCTION mol_in_surf2    
      
      CHARACTER(LEN=20) function str(k)
          !  "Convert an integer/real to string."
          REAL(KIND=8), INTENT(IN) :: k
          write (str, *) k
          str = adjustl(str)
      end function str
      
      FUNCTION nth(k, n )
          IMPLICIT NONE
          
          ! Declaring calling parameters:
          INTEGER, INTENT(IN) :: n            ! Length of string to return
          CHARACTER(len=20), INTENT(IN) :: k  ! String which is adjustl-ed
          CHARACTER(len=n) nth                ! Returned string
          ! Declaring local variables:
          
          ! Get string to return
          nth = k(1:n)
      END FUNCTION nth

      END MODULE 
MODULE traj_format2
!
! Purpose: 
! To declear data related to the simulation and traj.
IMPLICIT NONE
CONTAINS
  INTEGER FUNCTION sampling_number(nmo_start,nmo_end,ns)
    !
    ! Purpose:
    !  To calculate the total numbers of samples one want to include in their analysis.
    ! Data dictionary
    INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves

    write(*,*) 'In function sampling_number: nmo_end = ', nmo_end
    ! no. of samples = INT({no. of moves}/ns)
    positive: IF (nmo_end <0 .OR. nmo_start < 0 .OR. ns <0) THEN
      write(*,*) 'Please enter non-negative values for the ns, starting step and ending step.'
    ELSE IF (nmo_end < nmo_start) THEN
      write(*,*) 'Please note that starting step shoud not larger than  ending step.'
    ELSE IF (ns ==0) THEN
      sampling_number = nmo_end-(nmo_start-1)
    ELSE IF (nmo_end-(nmo_start-1) <= ns) THEN
      sampling_number = INT((nmo_end-(nmo_start-1))/ns + 1)
    ELSE IF (nmo_end-(nmo_start-1) > ns) THEN
      sampling_number = INT((nmo_end-(nmo_start-1))/ns)
    END IF positive
  END FUNCTION sampling_number

  SUBROUTINE read_traj_v1(indx,nmo_start,nmo_end,ns,nat,n_samples)
    ! Purpose:
    ! To read info from the trajectory file (format: ***.xyz)
    ! to READ data starting from a pattern-matched line.
    USE atom_module
    USE parameter_shared

    INTEGER :: iatom, imovie, i_sample
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    allocate(atom_info(nat,n_samples))
    i_sample = 1
    outer: DO imovie=1,nmo_end-nmo_start
      read(indx,*)!Neglect data of this line
      read(indx,120) sampled_movie(i_sample), sampled_time(i_sample)
      120 FORMAT (5X,I8,9X,F12.3,6X,F20.10)
      write(*,*) 'the step:', imovie
      inner: do iatom= 1,nat
        read (indx,*) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
          atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
        WRITE (*,*) & 
        atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), &
        atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
      enddo inner
      ! ns should be non-negative
      IF (ns < 0) THEN
        WRITE(*,*) 'Please note that ns should be non-negative.'
      ELSE IF (ns == 0) THEN
        CYCLE
      ELSE
        CALL skip_lines(indx, (nat+2)*(ns-1)) ! use CALL to run a SUBROUTINE
      ENDIF 

      i_sample = i_sample + 1

      ! To check if the sampling is finished
      check: IF (i_sample > n_samples) THEN 
        WRITE (*,*)'The total number of sample points are: ', n_samples
        EXIT
      END IF check
    ENDDO outer

  END SUBROUTINE read_traj_v1

  SUBROUTINE read_traj(indx,nmo_start,nmo_end,ns,nat,n_samples)
    ! Purpose:
    ! To read info from the trajectory file (format: ***.xyz)
    ! to READ data starting from a pattern-matched line.
    USE atom_module, ONLY: atom_info
    USE parameter_shared, ONLY: sampled_movie, sampled_time
    INTEGER :: iatom,i_sample
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER, INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    CHARACTER(LEN=4) :: head_char
    INTEGER :: y

    allocate(atom_info(nat,n_samples))
    i_sample = 1
    write(*,*) "read_traj(): New total time steps (n_samples):", n_samples
    DO WHILE (i_sample < n_samples+1) ! +1 means i_sample can take the value of n_samples 
        !read(indx, '(A4)') head_char
        read(indx, '(1X,A4)') head_char  ! for some other format, one can use this format
        PRE_CHECK:IF (head_char=="i = ") THEN
            BACKSPACE(UNIT=indx) ! Because I am not able to read other lines with the format '(A4,I8)', and have not find any good way, so I try to read it in '(A4)' first 
            !read(indx, '(A4,I8)') head_char, y
            read(indx, '(1X,A4,I8)') head_char, y
            CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-nmo_start,ns) == 0) THEN
            !Jie: For special case of ns=1, MOD(y-(nmo_start-1),ns) is always 0. Hence, it needs to be checked separately. 
            !Use ‘&’ to continue the line to avoid Fortran maximum line length of 132 characters. (Stupid Fortran!)
            !CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. &
            !                (ns == 1 .or. MOD(y-(nmo_start-1),ns) == 1))  THEN 
                BACKSPACE(UNIT=indx) ! Because we have to read the whole line with ' i = ' line.
                read(indx,130) sampled_movie(i_sample), sampled_time(i_sample)
                130 FORMAT (5X,I8,9X,F12.3,6X,F20.10)
                131 FORMAT (A4,3F20.10)
                inner: do iatom= 1,nat
                  read (indx,*) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
                    atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
                  if (atom_info(iatom, i_sample)%atom_name == "O") THEN
                      atom_info(iatom, i_sample)%mass = 16.00
                  elseif (atom_info(iatom, i_sample)%atom_name == "H") THEN
                      atom_info(iatom, i_sample)%mass = 1.00
                  endif
                enddo inner
                i_sample = i_sample + 1 !The position is important. It must be located before ENDIF MATCH
            ENDIF CHECK_HEAD
        ENDIF PRE_CHECK
    END DO
  END SUBROUTINE read_traj

  SUBROUTINE skip_lines(indx, i_input)
    !
    ! Purpose: 
    ! To skip lines when read data from the input
    IMPLICIT NONE
    INTEGER :: i
    INTEGER,INTENT(IN) :: i_input,indx
    do i=1,i_input
       read(indx,*) !Neglect (nat+2)*(ns-1) lines
    enddo    
  END SUBROUTINE skip_lines

END MODULE traj_format2
  MODULE types
  ! Derived data type to stroe array
  TYPE :: real_array
    CHARACTER(2) :: atom_name
    INTEGER :: time_step
    REAL(8), DIMENSION(3) :: coord_array
    TYPE (real_array), POINTER :: p
  END TYPE real_array
  ! Derived data type to store water molecules around a center atom
  TYPE :: water_array
    INTEGER :: time_step
    CHARACTER(2) :: center_atom_name
    REAL(8), DIMENSION(9) :: coord_water_array
    REAL(8) :: d_RO
    REAL(8), DIMENSION(2) :: theta 
    TYPE (water_array), POINTER :: p
  END TYPE water_array
  END MODULE types
MODULE water_module
! Purpose:
! To define the derived data type for water
USE atom_module

IMPLICIT NONE
TYPE :: water
  INTEGER :: molecule_id
  TYPE (atom):: oxygen
  TYPE (atom):: hydrogen1
  TYPE (atom):: hydrogen2
  REAL :: mass
END TYPE water

! The array water_info can be shared by subroutines  
TYPE(water), ALLOCATABLE, DIMENSION(:,:) :: water_info
END MODULE water_module
MODULE water_pair_module
!
! Purpose:
! To define the derived data type for a water pair

USE water_module
IMPLICIT NONE
TYPE :: water_pair
  INTEGER :: water_pair_id
  INTEGER :: head_id ! For example, 2
  INTEGER :: tail_id ! For example, 1
  TYPE (water) :: head ! water molecule
  TYPE (water) :: tail  !water molecule
  LOGICAL :: h_state 
  REAL :: h ! Hydrogen bond population operator 1 or 0 
  REAL :: hb_length 
END TYPE water_pair

!The array ```water_pair_info``` can be shared by subroutines  
TYPE(water_pair), ALLOCATABLE, DIMENSION(:,:) :: water_pair_info
END MODULE water_pair_module  
