
  REAL(KIND=8) FUNCTION distance2(r1,r2,boxsize)
      ! Date: 2021-3-24
      ! This is the correct version for calculate distance between two particles with PBC considered.
      USE module_ihb, ONLY: rk
      IMPLICIT NONE

      real(kind=rk), DIMENSION(3), INTENT(IN) :: r1,r2
      real(kind=rk), DIMENSION(3), INTENT(IN) :: boxsize
      REAL(kind=rk) :: dx,dy,dz
      REAL(kind=rk) :: temp

      dx =0; dy=0; dz=0; temp=0
 
      dx = r1(1) - r2(1)
      temp = boxsize(1)*0.5d0
      if (dx > temp) then
          dx = boxsize(1) - dx
      else if (dx < - temp) then
          dx = boxsize(1) + dx
      endif
      dy = r1(2) - r2(2)
      temp = boxsize(2)*0.5d0
      if (dy > temp) then
          dy = boxsize(2) - dy
      else if (dy < - temp) then
          dy = boxsize(2) + dy
      endif
      dz = r1(3) - r2(3)
      temp = boxsize(3)*0.5d0
      if (dz > temp) then
          dz = boxsize(3) - dz
      else if (dz < - temp) then
          dz = boxsize(3) + dz
      endif
      distance2 = dx**2 + dy**2 + dz**2
  END FUNCTION distance2
 
  REAL(kind=8) FUNCTION P2(x)
      IMPLICIT NONE
      ! P2: 2-order Legendre polynomial
      integer,parameter :: rk=8  
      real(kind=rk),INTENT(IN) :: x
      P2 = 0.5d0 * (3.d0 * x**2 -1.d0)
  END FUNCTION P2 

  REAL(kind=8) FUNCTION diff_axis (u1, u2, h)
      IMPLICIT NONE
      ! Considering PBC
      ! u2 is used as origin
      integer,parameter :: rk=8  
      real(kind=rk),INTENT(IN) :: u1,u2
      real(kind=rk),INTENT(IN) :: h
      real(kind=rk) :: du 
      du = u1-u2
      if (abs(du) > 0.5d0*h) then
          diff_axis = h - du 
      else
          diff_axis= du
      endif
  END FUNCTION diff_axis
  
  REAL(kind=8) FUNCTION pm_adh(r1,r2,r3,boxsize)
      IMPLICIT NONE
      INTEGER, PARAMETER :: rk=8  
      real(kind=rk), DIMENSION(3), INTENT(IN) :: r1,r2,r3
      real(kind=rk), DIMENSION(3), INTENT(IN) :: boxsize 
      REAL(kind=rk) :: diff_axis
      pm_adh=diff_axis(r3(1),r2(1),boxsize(1))*     &
             diff_axis(r1(1),r2(1),boxsize(1))+     &
             diff_axis(r3(2),r2(2),boxsize(2))*     &
             diff_axis(r1(2),r2(2),boxsize(2))+     &
             diff_axis(r3(3),r2(3),boxsize(3))*     &
             diff_axis(r1(3),r2(3),boxsize(3))      ! pm: point multiplication. 
  END FUNCTION pm_adh

  REAL(kind=8) FUNCTION pm_ahd(r1,r2,r3,boxsize)
      IMPLICIT NONE
      INTEGER, PARAMETER :: rk=8  
      REAL(kind=rk), DIMENSION(3), INTENT(IN) :: r1,r2,r3
      REAL(kind=rk), DIMENSION(3), INTENT(IN) :: boxsize 
      REAL(kind=rk) :: diff_axis
      pm_ahd=diff_axis(r1(1),r3(1),boxsize(1))*     &
             diff_axis(r2(1),r3(1),boxsize(1))+     &
             diff_axis(r1(2),r3(2),boxsize(2))*     &
             diff_axis(r2(2),r3(2),boxsize(2))+     &
             diff_axis(r1(3),r3(3),boxsize(3))*     &
             diff_axis(r2(3),r3(3),boxsize(3))      ! pm: point multiplication. 
  END FUNCTION pm_ahd 
  
  INTEGER FUNCTION get_total_number_of_lines(file_)
      IMPLICIT NONE
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

      get_total_number_of_lines = nlines
  END FUNCTION 

  FUNCTION get_ndx_of_oxygen(pos_filename, natoms) 
      implicit none
      !Parameters 
      INTEGER, PARAMETER :: rk=8

      !variables
      !=========
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER, DIMENSION(natoms) :: get_ndx_of_oxygen 
      !Local 
      integer :: i,nmovie,iatom,imovie
      real(kind=rk), allocatable, dimension (:,:)  :: x,y,z
      character(LEN=3), allocatable, dimension (:) :: atom_type
      
      !Initialization
      get_ndx_of_oxygen = 0
      !==================
      !Read data in input
      !==================
      nmovie=1 ! Number of movie steps reqired.
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
                  get_ndx_of_oxygen(i) = iatom 
              endif
          enddo
      enddo
      close(10)
      deallocate(atom_type,x,y,z)
  END FUNCTION get_ndx_of_oxygen

  FUNCTION get_number_of_oxygen(pos_filename, natoms) 
      implicit none
      !Parameters 
      INTEGER, PARAMETER :: rk=8

      !variables
      !=========
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER :: get_number_of_oxygen 
      !Local 
      integer :: i,nmovie,iatom,imovie
      real(kind=rk), allocatable, dimension (:,:)           :: x,y,z
      character(LEN=3), allocatable, dimension (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1 ! Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  ! Note that nmovie=1
          read(10,*) ! Neglect data of this line
          read(10,*)                  
          i=0
          do iatom= 1, natoms
              read (10,*) atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O') then
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_oxygen = i
      close(10)
      deallocate(atom_type,x,y,z)
  END FUNCTION get_number_of_oxygen

  FUNCTION get_number_of_nitrogen(pos_filename, natoms) 
      implicit none
      !Parameters and variables
      !========================
      integer,parameter :: rk=8  
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER :: get_number_of_nitrogen 
      !Local 
      integer :: i, nmovie, iatom, imovie
      real(kind=rk), allocatable, dimension (:,:)           :: x, y, z
      character(LEN=3), allocatable, dimension (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1 ! Number of movie steps reqired.
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
              if (trim(atom_type(iatom)) .eq. 'N') then
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_nitrogen = i
      close(10)
      deallocate(atom_type,x,y,z)
  END FUNCTION get_number_of_nitrogen

  FUNCTION get_number_of_lithium(pos_filename, natoms) 
      implicit none
      !Parameters and variables
      !========================
      integer,parameter :: rk=8  
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER :: get_number_of_lithium 
      !Local 
      integer :: i, nmovie, iatom, imovie
      real(kind=rk), allocatable, dimension (:,:)           :: x, y, z
      character(LEN=3), allocatable, dimension (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1 ! Number of movie steps reqired.
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
              if (trim(atom_type(iatom)) .eq. 'Li') then
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_lithium = i
      close(10)
      deallocate(atom_type,x,y,z)
  END FUNCTION get_number_of_lithium

  FUNCTION hydrogen_ndx_list(ndx_oxygen_1, & 
           ndx_oxygen_2,pos_filename,natoms,boxsize) 
      implicit none
      ! 2020/07
      !========================
      !Parameters and variables
      !========================
      integer,parameter :: rk=8  
      INTEGER,INTENT(IN) :: ndx_oxygen_1, ndx_oxygen_2
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      real(kind=rk),dimension(3),INTENT(IN) :: boxsize
      
      !Local 
      REAL(kind=rk) :: distance2
      INTEGER :: get_number_of_oxygen,get_number_of_hydrogen
      integer :: i,j,iatom,nmovie,& 
                    m1,m2,m3,i_H,&
                    i1,i2,ii,jj,i_OW, num
      INTEGER, DIMENSION(4) :: hydrogen_ndx_list
      real(kind=rk), allocatable,dimension (:,:) :: x,y,z
      character(LEN=3), allocatable, dimension (:) :: atom_type
      integer, allocatable, dimension (:) :: ndx_OW,ndx_H
      real(kind=rk), dimension(3) :: r_a, r_b
      real(kind=rk), parameter :: r_ohc=1.21d0   ! rOH (1.1**2)
      real(kind=rk) :: r ! distance between O and H.
      ! Initialization
      i1=0; i2=0
      !==================
      !Read data in input
      !==================
      nmovie=1 !Number of movie steps reqired.
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
              r_a= (/x(m2,j),y(m2,j),z(m2,j)/) ! Coordinate of O2
              r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
              r=distance2(r_a,r_b,boxsize) ! squared distance between O2 and H
              if (r<r_ohc) then ! If the H is bound to O2
                  num = num + 1
                  hydrogen_ndx_list(num) = m3 ! Add the index of the H atom into the H list for the water pair (O1, O2)
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
      integer,parameter :: rk=8              
      INTEGER, INTENT(IN) :: ndx_X, ndx_oxygen_2
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      real(kind=rk), dimension(3), INTENT(IN) :: boxsize
      
      !Local 
      real(kind=rk), dimension(3) :: r_a, r_b
      real(kind=rk),parameter :: r_ohc=1.21d0  ! rOH (1.1**2)
      real(kind=rk) :: r ! distance between O and H.
      integer :: i,j,nmovie,iatom,& 
              m1,m2,m3,i_H,&
              i1,i2,ii,jj,i_O,num
      INTEGER, DIMENSION(2) :: hydrogen_ndx_list_XO ! the length is 2, because there are 2 Hydrogen are bonded to ndx_oxygen_2 but no Hydrogen is bonded to ndx_X
      real(kind=rk),allocatable,dimension (:,:) :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:) :: ndx_O,ndx_H
      real(kind=rk) :: distance2
      integer :: get_number_of_hydrogen, get_number_of_oxygen
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
      !=======================================================
      ! Calculate the indices of O (H) atoms in water molecules
      !=======================================================
      allocate(ndx_O(i_O)) ! this should be put after i_O is defined
      allocate(ndx_H(i_H)) ! this should be put after i_H is defined
      !allocate(ndx_I(i_I)) ! this should be put after i_I is defined
      i=0; ii=0; jj=0
      !write(*,*) "pos_filenae:",pos_filename
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
                 !write(*,*) "iatom:(test O) ", iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ndx_H(ii+1)=iatom
                 ii=ii+1
                 !write(*,*) "iatom:(test) ", iatom
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
              !write(*,*) "i_H: ", i_H
              !write(*,*) "m3: ", m3
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
  END FUNCTION hydrogen_ndx_list_XO

  FUNCTION get_number_of_hydrogen(pos_filename, natoms) 
      implicit none
      !========================
      integer,parameter :: rk=8              
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      INTEGER :: get_number_of_hydrogen
      !Local 
      integer :: i,nmovie,iatom,imovie
      real(kind=rk),allocatable,dimension (:,:) :: x,y,z
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
  END FUNCTION get_number_of_hydrogen

  FUNCTION get_number_of_iodine(pos_filename, natoms) 
      implicit none
      integer,parameter :: rk=8              
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      INTEGER :: get_number_of_iodine 
      !Local 
      integer :: i, nmovie, iatom, imovie
      real(KIND=rk),allocatable,dimension (:,:) :: x,y,z
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
  END FUNCTION get_number_of_iodine

  !INTEGER FUNCTION grid_index(x,y,divx,divy,nb_divx,nb_divy)
  !    ! transfer the coordinates (x,y) to grid_index, which is an integer
  !    implicit none
  !    integer,parameter :: rk=8              
  !    REAL(kind=rk), INTENT(IN) :: x,y
  !    REAL(kind=rk), INTENT(IN) :: divx, divy
  !    INTEGER, INTENT(IN) :: nb_divx, nb_divy
  !    INTEGER, DIMENSION(2) :: ind
  !    
  !    !Initialization
  !    ind = 0

  !    ind(1) = floor(x/divx) 
  !    ind(2) = floor(y/divy)

  !    !Correction for avoiding 0 index
  !    !IF (ind(1) == 0) THEN
  !    !    ind(1) = ind(1) + 1
  !    !ENDIF
  !    !IF (ind(2) == 0) THEN
  !    !    ind(2) = ind(2) + 1
  !    !ENDIF
  !    if (ind(1) == 0) then
  !      grid_index = ind(2) 
  !    else
  !      grid_index = ind(1)*nb_divy + ind(2) 
  !    endif
  !    grid_index = grid_index + 1
  !END FUNCTION grid_index

  INTEGER FUNCTION grid_index(x,y,divx,divy,n_divy)
      ! transfer the coordinates (x,y) to grid_index, which is an integer
      INTEGER, PARAMETER :: rk=8
      REAL(KIND=rk), INTENT(IN) :: x,y
      REAL(KIND=rk), INTENT(IN) :: divx, divy
      INTEGER, INTENT(IN) :: n_divy
      INTEGER, DIMENSION(2) :: ind

      !Initialization
      ind = 0.0

      ind(1) = FLOOR(x/divx)
      ind(2) = FLOOR(y/divy)

      grid_index = ind(1) * n_divy + ind(2)+1

  END FUNCTION grid_index


  LOGICAL FUNCTION oh_in_surf1(surf1_mol1,z1,thickness)
      implicit none
      INTEGER, PARAMETER :: rk=8
      LOGICAL :: mol1_in_surf1
      REAL(kind=rk), INTENT(IN) :: surf1_mol1
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z1
      mol1_in_surf1 = surf1_mol1 + & 
          thickness > z1
      oh_in_surf1 = mol1_in_surf1 
  END FUNCTION oh_in_surf1    

  LOGICAL FUNCTION oh_in_surf2(surf2_mol1,z1,thickness)
      implicit none
      INTEGER, PARAMETER :: rk=8
      LOGICAL :: mol1_in_surf2
      REAL(kind=rk), INTENT(IN) :: surf2_mol1
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z1
      mol1_in_surf2 = surf2_mol1 - & 
          thickness < z1
      oh_in_surf2 = mol1_in_surf2 
  END FUNCTION oh_in_surf2    

  LOGICAL FUNCTION pair_in_surf1(surf1_mol1,z1,surf1_mol2,z2,thickness)
      implicit none
      INTEGER, PARAMETER :: rk=8
      LOGICAL :: mol1_in_surf1, mol2_in_surf1
      REAL(kind=rk), INTENT(IN) :: surf1_mol1,surf1_mol2
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z1,z2
      mol1_in_surf1 = surf1_mol1 + & 
          thickness > z1
      mol2_in_surf1 = surf1_mol2 + &
          thickness > z2
      pair_in_surf1 = mol1_in_surf1 .and. mol2_in_surf1
  END FUNCTION pair_in_surf1    

  LOGICAL FUNCTION pair_in_surf2(surf2_mol1,z1,surf2_mol2,z2,thickness)
      implicit none
      INTEGER, PARAMETER :: rk=8
      LOGICAL :: mol1_in_surf2, mol2_in_surf2
      REAL(kind=rk), INTENT(IN) :: surf2_mol1,surf2_mol2
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z1,z2
      mol1_in_surf2 = surf2_mol1 - & 
          thickness < z1
      mol2_in_surf2 = surf2_mol2 - &
          thickness < z2
      pair_in_surf2 = mol1_in_surf2 .and. mol2_in_surf2
  END FUNCTION pair_in_surf2    

  LOGICAL FUNCTION atom_in_surf1(surf1_mol1,z1,thickness)
      INTEGER, PARAMETER :: rk=8
      REAL(KIND=rk), INTENT(IN) :: surf1_mol1
      REAL(KIND=rk), INTENT(IN) :: z1,thickness
      LOGICAL :: mol1_in_surf1

      mol1_in_surf1 = surf1_mol1 + &
          thickness > z1
      atom_in_surf1 = mol1_in_surf1

  END FUNCTION atom_in_surf1

  LOGICAL FUNCTION atom_in_surf2(surf2_mol1,z1,thickness)
      INTEGER, PARAMETER :: rk=8    
      REAL(KIND=rk), INTENT(IN) :: surf2_mol1
      REAL(KIND=rk), INTENT(IN) :: z1,thickness
      LOGICAL :: mol1_in_surf2

      mol1_in_surf2 = surf2_mol1 - &
          thickness < z1
      atom_in_surf2 = mol1_in_surf2

  END FUNCTION atom_in_surf2

  LOGICAL FUNCTION mol_in_surf1(surf1_mol,z1,thickness)
      !Check if a molecule is in the surf1 ( the lower layer)
      implicit none
      INTEGER, PARAMETER :: rk=8
      REAL(kind=rk), INTENT(IN) :: surf1_mol
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z1
      mol_in_surf1 = surf1_mol + thickness > z1
  END FUNCTION mol_in_surf1    

  LOGICAL FUNCTION mol_in_surf2(surf2_mol,z2,thickness)
      !Check if a molecule is in the surf2 ( the upper layer)
      implicit none
      INTEGER, PARAMETER :: rk=8
      REAL(kind=rk), INTENT(IN) :: surf2_mol
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z2
      mol_in_surf2 = surf2_mol - thickness < z2
  END FUNCTION mol_in_surf2    
  
  character(len=20) function int2str(k)
      !  "Convert an integer to string."
      integer, intent(in) :: k
      write (int2str, *) k
      int2str = adjustl(int2str)
  end function int2str

  character(len=20) function str(k)
      !  "Convert an integer/real to string."
      real (kind=8), intent(in) :: k
      write (str, '(F3.1)') k
      str = adjustl(str)
  end function str

  FUNCTION nth(k,n)
      ! To return a string containing the first N characters
      ! of the alphabet.
      IMPLICIT NONE
  
      ! Declaring calling parameters:
      INTEGER, INTENT(IN) :: n ! Length of string to return
      CHARACTER(len=20), INTENT(IN) :: k ! String which is adjustl-ed
      CHARACTER(len=n) nth ! Returned string
      
      ! Get string to return
      nth = k(1:n)
  END FUNCTION nth

  INTEGER FUNCTION sampling_number(nmo_start,nmo_end,ns)
      ! To calculate the total numbers of samples one WANT to include in their analysis.
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
      INTEGER, INTENT(IN) :: ns ! Get one sample from the trajectory every ns step.

      positive: IF (nmo_end <0 .OR. nmo_start < 0 .OR. ns <0) THEN
        write(*,*) 'Please enter non-negative values for the ns, starting step and ending step.'
      ELSE IF (nmo_end < nmo_start) THEN
        write(*,*) 'Please note that starting step shoud not larger than  ending step.'
      ELSE IF (ns == 0) THEN
        sampling_number = nmo_end-(nmo_start-1)
      ELSE 
        sampling_number = FLOOR(FLOAT(nmo_end-(nmo_start-1))/FLOAT(ns))
      END IF positive
  END FUNCTION sampling_number

  SUBROUTINE read_traj_3(indx,nmo_start,nmo_end,ns,nat,n_samples,sampled_movie,sampled_time,boxsize,criterion)
      ! Purpose:
      !  To read info from the trajectory file (format: ***.xyz)
      !  1. Read data starting from a pattern-matched line.
      !  2. Determine free OH groups.
      USE atom_module, ONLY: atom, atom_info
      USE water_molecule_types, ONLY: hydrogen_atom, oxygen_atom, O_info, H_info
      IMPLICIT NONE

      REAL,PARAMETER :: rooc=12.25 ! cutoff distance of rOO (3.5**2 )
      REAL,PARAMETER :: cosPhiC123=0.866 ! 1.732/2; phiC=pi/6.
      REAL,PARAMETER :: cosPhiC123_freeOH= 0.6428 !  phiC= 50 degree. (Ref. J. Chem. Theory Comput. 2018, 14, 357−364)
      REAL,PARAMETER :: cosPhiC132=-0.5 ! -1./2; phiC132=2pi/3.
      REAL,PARAMETER :: cosPhiC132_freeOH = -0.342 ! ; phiC132= 110. (Ref. Tang, J. Chem. Theory Comput. 2018, 14, 357−364)
      INTEGER, PARAMETER :: rk=8
      REAL(KIND=rk) :: distance2, pm_adh, pm_ahd
      !REAL(kind=rk) :: distance2, pm_adh, pm_ahd
      REAL(KIND=rk),PARAMETER :: max_time_for_corr = 12.0 ! Unit: ps. 
      REAL(KIND=rk),PARAMETER :: h_min=0.5 ! condition for the existence of a h-bond for a step
      REAL(KIND=rk),PARAMETER :: hb_min=0.5 ! condition for the existence of h-bond for a pair of water molecules
      REAL(KIND=rk),PARAMETER :: nfb_min=0.5 ! condition for the existence of free OH for a  OH in water molecule

      CHARACTER(LEN=4) :: head_char
      INTEGER :: iatom, i_sample
      INTEGER :: i_O, i_H
      INTEGER :: n_H ! number of OH groups; or num of H atoms
      INTEGER :: n_O ! number of O atoms
      INTEGER :: i, j, k, jj, bond
      INTEGER :: k_H, k_O, k_O1, k_O2
      INTEGER,INTENT(IN) :: criterion 
      INTEGER,INTENT(IN) :: indx
      INTEGER,INTENT(IN) :: nat
      INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
      INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
      INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
      REAL(KIND=rk),DIMENSION(3),INTENT(IN) :: boxsize      
      REAL(KIND=rk) :: r13, cosphi, pm, cosphi_, pm_, norm_rr
      REAL(KIND=rk) :: r21, r31, r32, r23 ! For the second criterion of HB

      REAL(KIND=rk),DIMENSION(3) :: r1, r2, r3 ! pbc 
      INTEGER :: m1, m2, m3
      INTEGER :: idx_O1 ! The self index of O1 in all O atoms
      INTEGER :: idx_H ! The self Index of H in all H atoms

      !TYPE(atom),DIMENSION(nat,n_samples),INTENT(INOUT) :: atom_info
      INTEGER,DIMENSION(n_samples) :: sampled_movie
      REAL(kind=rk),DIMENSION(n_samples) :: sampled_time
      INTEGER :: y
      INTEGER :: tmp_index
      INTEGER :: num_limit ! Temperary varible, represents the largest number of an O atoms' neighbors.
      LOGICAL,ALLOCATABLE,DIMENSION (:) :: is_free

      !allocate(atom_info(nat,n_samples)) ! Define atom_info

      n_O = nat / 3
      n_H = nat * 2 / 3
      i_sample = 1
      DO WHILE (i_sample < n_samples+1) ! +1 means i_sample can take the value of n_samples 
          read(indx, '(1X,A4)') head_char  ! for some other format, one can use this format
          PRE_CHECK:IF (head_char=="i = ") THEN
              BACKSPACE(UNIT=indx) ! Because I am not able to read other lines with the format '(A4,I8)', and have not find any good way, so I try to read it in '(A4)' first 
              read(indx, '(1X,A4,I8)') head_char, y ! for some other format, one can use this format
              CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-nmo_start,ns) == 0) THEN
                  !-------------------------------------------------------------------------------------------------------
                  !We use y>nmo_start-1, because we want to consider the first step 'i=0'
                  !-------------------------------------------------------------------------------------------------------
                  BACKSPACE(UNIT=indx) ! Because we have to read the whole line with ' i = ' line.
                  !read(indx,130) sampled_movie(i_sample), sampled_time(i_sample), sampled_energy(i_sample)
                  read(indx,120) sampled_movie(i_sample), sampled_time(i_sample)
                  130 FORMAT (1X,4X,I8,9X,F12.3,6X,F20.10)
                  120 FORMAT (1X,4X,I8,9X,F12.3)
                  131 FORMAT (A4,3F20.10)
                  i_O = 0
                  i_H = 0
                  inner: do iatom= 1,nat
                    read (indx,*) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
                      atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
                    if (atom_info(iatom, i_sample)%atom_name == "O") THEN
                        i_O = i_O + 1
                        atom_info(iatom, i_sample)%mass = 16.00
                        atom_info(iatom, i_O)%atom_name = "O"
                        O_info(i_O,i_sample)%coord = atom_info(iatom,i_sample)%coord
                        O_info(i_O,i_sample)%atom_id = iatom
                        O_info(i_O,i_sample)%host_id = iatom
                        O_info(i_O,i_sample)%H_ids(1) = iatom + 1
                        O_info(i_O,i_sample)%H_ids(2) = iatom + 2
                        O_info(i_O,i_sample)%num_oxygen_neighbors = 0
                        O_info(i_O,i_sample)%self_indices_oxygen_neighbors = 0
                        O_info(i_O,i_sample)%indices_oxygen_neighbors = 0
                    elseif (atom_info(iatom, i_sample)%atom_name == "H") THEN
                        i_H = i_H + 1
                        atom_info(iatom, i_sample)%mass = 1.00
                        atom_info(i_H, i_sample)%atom_name = "H"
                        H_info(i_H, i_sample)%coord = atom_info(iatom,i_sample)%coord
                        H_info(i_H, i_sample)%atom_id = iatom
                        ! Assgining host_id. 3: water molecule has 3 atoms.
                        if ( mod(iatom,3) == 2 ) then
                            H_info(i_H, i_sample)%host_id = iatom - 1
                        else if ( mod(iatom,3) == 0 ) then
                            H_info(i_H, i_sample)%host_id = iatom - 2
                        else
                           write(*,*) "The order of O and H is not `O H H O H H...`."
                        endif
                        H_info(i_H, i_sample)%IsFree = .False.
                    endif
                  enddo inner
                  i_sample = i_sample + 1 !The position is important. It must be located before ENDIF 
              ENDIF CHECK_HEAD
          ENDIF PRE_CHECK
      END DO
      write(*,*) "Total numFrames", i_sample -1 ! To check the total number of samples)


      !A loop over TOTAL Time Step. Purpose: Find Oxygen atoms' Oxygen neighbors.
      TTLOOP: DO jj = 1,n_samples
      DO k_O1 = 1, n_O
          ! For a given O, we calculate the list of Oxygen neighbors of this O atom. 
          tmp_index = 0
          m1 = 3 * k_O1 - 2 ! Host O1, total index
          !write(*,*) "m1=", m1
          DO k_O2 = 1 , n_O
              if (k_O2 .ne. k_O1) then
                  m2 = 3 * k_O2 - 2 ! Total order of O2. For searching neighbors of Host O1
                  ! The total inices of H atoms related to O1 (Host O) and O2 (potential neighor O)
                  ! Check if O2 is a neighbor of O1
                  r1 = (/atom_info(m1,jj)%coord(1),atom_info(m1,jj)%coord(2),&
                     atom_info(m1,jj)%coord(3) /)
                  r2 = (/atom_info(m2,jj)%coord(1),atom_info(m2,jj)%coord(2),&
                         atom_info(m2,jj)%coord(3) /)
                  r21 = distance2(r1, r2, boxsize)
                  if ((r21 .lt. rooc)) THEN
                     write(*,*) "m2=", m2
                     ! In the neighbor list of O1, add the total indices of O2.
                     O_info(k_O1,jj)%num_oxygen_neighbors = O_info(k_O1,jj)%num_oxygen_neighbors + 1 ! increase the number of O neighbors
                     tmp_index = O_info(k_O1,jj)%num_oxygen_neighbors
                     write(*,*) "num_oxygen_neighbors for",k_O1,"-th O:", jj, "-th step", tmp_index
                     O_info(k_O1,jj)%self_indices_oxygen_neighbors(tmp_index) = k_O2
                     O_info(k_O1,jj)%indices_oxygen_neighbors(tmp_index) = m2
                     O_info(k_O2, jj )%atom_id = m2 ! define total index of O2
                  endif
              endif
          ENDDO
      ENDDO
      ENDDO TTLOOP


      !The purpose of this loop is to  define the free OH groups at each time step.
      tLOOP: DO jj = 1,n_samples
          !===================================
          !The loop to find out free OH groups
          !===================================
          Ohost: DO k_O1 = 1, n_O ! self index of O1
              idx_O1 = k_O1 * 3 - 2 ! The total index of O1
              m1 = idx_O1 ! Indices of total index of the Host 
              OH: DO bond = 1,2
                  idx_H = O_info(k_O1,1)%H_ids(bond) ! 1: the first time step
                  m3 = idx_H 
                  if (bond == 1) then
                      k_H = (idx_H + 1) * 2/3 - 1
                  else
                      k_H = (idx_H ) * 2/3 
                  endif
                  if (m1 < 1) then
                      write(*,*) "Error: m1 must be positive!"
                  endif
                  if (O_info(k_O1,jj)%num_oxygen_neighbors .le. 0) then
                      H_info(k_H,jj)%IsFree = .True.
                  else
                      num_limit = O_info(k_O1,jj)%num_oxygen_neighbors
                      allocate(is_free(num_limit))
                      DO k_O2 = 1, num_limit 
                          m2 = O_info(k_O1,jj)%indices_oxygen_neighbors(k_O2) ! Indices of total index of the neighbor O2 
                          ! For all tuples (idx_O1, idx_O_neighbor[1]), (idx_O1, idx_O_neighbor[2]), (idx_O1, idx_O_neighbor[3]),...
                          !Check whether the OH is free OH group respectively. Then based on these results, we can know that whether
                          !the OH is free OH or not, so that can Calculate nf(j) later.

                          r1 = (/atom_info(m1,jj)%coord(1),atom_info(m1,jj)%coord(2),&
                                 atom_info(m1,jj)%coord(3) /)
                          r2 = (/atom_info(m2,jj)%coord(1),atom_info(m2,jj)%coord(2),&
                                 atom_info(m2,jj)%coord(3) /)
                          r3 = (/atom_info(m3,jj)%coord(1),atom_info(m3,jj)%coord(2),&
                                 atom_info(m3,jj)%coord(3) /)
                          r21 = distance2(r1, r2, boxsize) 
                          if (criterion == 1) THEN
                              r23 = distance2(r3, r2, boxsize) 
                              r13 = distance2(r3, r1, boxsize) 
                              pm_ = pm_adh(r2,r1,r3,boxsize) ! if H is bound to O1
                              !norm_rr = sqrt(r21*r23) ! if H is bound to O2 (Not for here.)
                              norm_rr = sqrt(r21*r13) ! if H is bounded to O1
                              cosphi_= pm_/norm_rr
                              if (cosphi_  .lt. cosPhiC123_freeOH)  THEN
                                  is_free(k_O2) = .True. ! The k_O2-th step for checking whether the k_H -th H is free.
                              endif
                          elseif (criterion == 2) THEN
                              !Follow the second cirterion of HB.
                              r31 = distance2(r1,r3,boxsize) 
                              r32 = distance2(r2,r3, boxsize) 
                              pm = pm_adh(r1,r2,r3,boxsize)
                              cosphi = pm/(sqrt(r31*r32))

                              !Follow the scond criterion of HB.
                              !-0.342 comes from cos(110 degree)
                              if (cosphi .gt. cosPhiC132_freeOH) THEN 
                                  is_free(k_O2) = .True.
                              endif
                          endif
                          if ( all(is_free .eqv. .True.) .eqv. .True. ) then
                              H_info(k_H,jj)%IsFree = .True.
                          endif
                      ENDDO
                      deallocate(is_free)
                  endif          
              ENDDO OH
          ENDDO Ohost
      ENDDO tLOOP
  END SUBROUTINE read_traj_3

  SUBROUTINE read_surf_coord(indx,n_samples,n_grid,surf_info_fortran)
    IMPLICIT NONE

    INTEGER :: ierror, i_sample, i_grid
    INTEGER :: skip_num ! The number of steps skipped before reading
    INTEGER,INTENT(IN) :: indx
    INTEGER,INTENT(IN) :: n_grid ! n_grid = nb_divx * nb_divy; KEEP
    INTEGER,INTENT(IN) :: n_samples ! n_samples = INT(nmo/ns); KEEP
    REAL(kind=8),DIMENSION(2,n_samples,n_grid),INTENT(INOUT) :: surf_info_fortran
    
    ierror =0; i_sample=0; i_grid=0; skip_num=0
  
    131 FORMAT (11X,2F13.6)
      outer: DO i_sample = 1, n_samples 
          read(indx, '(A4,I5)',IOSTAT=ierror) 
          inner: do i_grid = 1, n_grid
              read (indx,131,IOSTAT=ierror) surf_info_fortran(1,i_sample,i_grid), & 
                  surf_info_fortran(2,i_sample,i_grid)
          enddo inner
      END DO outer
  END SUBROUTINE read_surf_coord

  SUBROUTINE skip_lines(indx, i_input)
      ! To skip lines when read data from the input
      IMPLICIT NONE
      INTEGER :: i
      INTEGER,INTENT(IN) :: i_input,indx
      do i=1,i_input
         read(indx,*) !Neglect (nat+2)*(ns-1) lines
      enddo    
  END SUBROUTINE skip_lines

  FUNCTION hydrogen_ndx_list_for_oxygen(ndx_oxygen_1, & 
           pos_filename,natoms,boxsize) 
      implicit none
      ! 2020/07
      !========================
      !Parameters and variables
      !========================
      integer,parameter :: rk=8  
      INTEGER,INTENT(IN) :: ndx_oxygen_1
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      real(kind=rk),dimension(3),INTENT(IN) :: boxsize
      
      !Local 
      REAL(kind=rk) :: distance2
      INTEGER :: get_number_of_oxygen,get_number_of_hydrogen
      integer :: i,j,iatom,nmovie,& 
                    m1,m3,i_H,&
                    i1,ii,jj,i_OW, num
      INTEGER, DIMENSION(2) :: hydrogen_ndx_list_for_oxygen
      real(kind=rk), allocatable,dimension (:,:) :: x,y,z
      character(LEN=3), allocatable, dimension (:) :: atom_type
      integer, allocatable, dimension (:) :: ndx_OW,ndx_H
      real(kind=rk), dimension(3) :: r_a, r_b
      real(kind=rk), parameter :: r_ohc=1.21d0   ! rOH (1.1**2)
      real(kind=rk) :: r ! distance between O and H.
      ! Initialization
      i1=0
      !==================
      !Read data in input
      !==================
      nmovie=1 !Number of movie steps reqired.
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
      num = 0 ! For counting
      ! 
      do j =1, 1
          do ii=1,i_H
              m3=ndx_H(ii)
              ! I use pbc to consider the distance r
              r_a= (/x(m1,j),y(m1,j),z(m1,j)/) ! Coordinate of O
              r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
              r=distance2(r_a,r_b,boxsize) ! squared distance between O and H
              if (r<r_ohc) then ! If the H is bound to the O atom
                  num = num + 1
                  hydrogen_ndx_list_for_oxygen(num) = m3 ! Add the index of the H atom into the H list for the water pair (O1, O2)
              endif
          enddo
      enddo
      deallocate(ndx_OW,ndx_H,x,y,z)
      !=============================
  END FUNCTION hydrogen_ndx_list_for_oxygen
   
  SUBROUTINE make_neighbors(nat,n_samples,boxsize,neighbors_number,neighbors_indices,sphere_info)
    USE module_ihb, only: rk, sphere, num_coord_max, r_OH_LU2, r_ON_LU2, distance2
    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER,INTENT(IN) :: nat,n_samples
    real(kind=rk), DIMENSION(3), INTENT(IN) :: boxsize
    TYPE(sphere),DIMENSION(nat,n_samples),INTENT(INOUT) :: sphere_info
    INTEGER,DIMENSION(nat),INTENT(INOUT) :: neighbors_number
    INTEGER,DIMENSION(num_coord_max,nat),INTENT(INOUT) :: neighbors_indices
    
    INTEGER,PARAMETER :: start_step = 1
    REAL(kind=rk) :: r_LU2 ! LU: Least Upper bound; 2: squared value 

    !Initialize arrays and variables
    neighbors_number = 0
    neighbors_indices = 0
    r_LU2 = MAX(r_OH_LU2,r_ON_LU2)

    Outer: DO i = 1, nat
      Inner: DO j = 1, nat
        IF (j==i) CYCLE
        IF ( distance2(sphere_info(i,start_step)%coord, sphere_info(j,start_step)%coord, boxsize)<   &
           r_LU2                                                                                     &
           ) THEN 
           neighbors_number(i) = neighbors_number(i) + 1
           neighbors_indices(neighbors_number(i),i) = j
           !WRITE(*,*) "Atom Name:",sphere_info(i,start_step)%name
           !WRITE(*,*) "Atom id:",sphere_info(i,start_step)%id
           !WRITE(*,*) "distance:",SQRT(distance2(sphere_info(i,start_step)%coord, sphere_info(j,start_step)%coord, boxsize))
        END IF
      END DO Inner
    END DO Outer
  END SUBROUTINE make_neighbors

  SUBROUTINE obtain_indices(atom_name,indices,nat,n_samples,sphere_info)
    !===================================================================
    ! This subroutine does one thing for an chemical element (X, a string denoted as atom_name):
    ! obtain the indices of atom X and store it in the array indices
    ! indices: the output array
    !===========================
    
    USE module_ihb, only: sphere, num_coord_max
    IMPLICIT NONE
 
    CHARACTER(len=*),INTENT(IN) :: atom_name
    INTEGER,INTENT(IN) :: nat,n_samples
    TYPE(sphere),DIMENSION(nat,n_samples),INTENT(INOUT) :: sphere_info
    INTEGER,DIMENSION(nat), INTENT(INOUT) :: indices
    
    INTEGER,PARAMETER :: start_step=1 
    INTEGER :: i, n_count
    
    !WRITE(*,*) "obtain_indices:  atom_name = ", atom_name
    ! Initialize indices
    indices = 0
    ! Generate the array of index of Li,I.
    n_count = 0
    DO i = 1, nat
      IF ( trim(sphere_info(i,start_step)%name) == trim(atom_name) ) THEN
        n_count = n_count + 1
        indices(n_count) = i
        !WRITE(*,*)"indices(n_count): ", indices(n_count)
      ENDIF
    END DO
  END SUBROUTINE obtain_indices

  !SUBROUTINE make_wat_neighbors_for_X_shell(atom_name,step_t,nat,n_samples,boxsize, &
  !           wat_neighbors_number,wat_neighbors_indices,sphere_info)
  !  !===================================================================
  !  ! this subroutine do one thing:
  !  ! For an atom (index = i_ion), at time t (denoted by step_t), 
  !  ! 1: to find the number (N) of its neighboring water molecules and store N in wat_neighbors_number(i_ion); 
  !  ! 2: to add the index of O of this water molecule into the wat_neighbors_indices(ind, i_ion)
  !  !====================================================================
  !  
  !  USE module_ihb, only: rk,sphere,num_coord_mol_max,r_OH_LU2,r_NN_OW_LU2,r_ON_OW_LU2,r_Li_OW_LU2,&
  !    r_Na_OW_LU2,r_K_OW_LU2,r_I_OW_LU2,r_OW_OW_LU2, distance2
  !  IMPLICIT NONE

  !  INTEGER :: i, j, k
  !  !INTEGER :: n_count  ! To count the number of neighbors of X-ion/atom/etc.
  !  CHARACTER(len=*),INTENT(IN) :: atom_name
  !  INTEGER,INTENT(IN) :: nat,n_samples
  !  real(kind=rk), DIMENSION(3), INTENT(IN) :: boxsize
  !  TYPE(sphere),DIMENSION(nat,n_samples),INTENT(INOUT) :: sphere_info
  !  INTEGER,DIMENSION(nat),INTENT(INOUT) :: wat_neighbors_number
  !  INTEGER,DIMENSION(num_coord_mol_max,nat),INTENT(INOUT) :: wat_neighbors_indices
  !   
  !  ! Local variable
  !  INTEGER,DIMENSION(nat) :: indices
  !  !REAL(kind=rk) :: dist2
  !  
  !  INTEGER, INTENT(IN):: step_t 
  !  REAL(kind=rk) :: r_LU2 ! LU: Least Upper bound; 2: squared value 

  !  !Initialize arrays and variables
  !  i=0;j=0;k=0
  !  wat_neighbors_number = 0
  !  wat_neighbors_indices = 0

  !  !================    
  !  ! Determine r_LU2
  !  !================    
  !  IF (atom_name =="Li") THEN
  !    r_LU2 = r_Li_OW_LU2
  !  ELSE IF (atom_name =="Na") THEN
  !    r_LU2 = r_Na_OW_LU2
  !  ELSE IF (atom_name =="K") THEN
  !    r_LU2 = r_K_OW_LU2
  !  ELSE IF (atom_name =="I") THEN
  !    r_LU2 = r_I_OW_LU2
  !  ELSE IF (atom_name =="N") THEN
  !    r_LU2 = r_NN_OW_LU2
  !  END IF

  !  ! Generate the array of index of Li, or I, etc.
  !  WRITE(*,*) "make_wat_neighbors_for_X_shell:  atom_name = ", atom_name
  !  CALL obtain_indices(atom_name,indices,nat,n_samples,sphere_info)
  !  !Test
  !  !WRITE(*,*) "indices of atom_name: ", indices
  !  Outer: DO i = 1, nat
  !    IF(indices(i) == 0) EXIT
  !    k = indices(i)
  !    WRITE(*,*) "i=",i, "; indices(i)=", k
  !    Inner: DO j = 1, nat
  !    IF (j== k ) CYCLE
  !    IF (TRIM(sphere_info(j,step_t)%name)=="O" .AND. &
  !        distance2(sphere_info(k,step_t)%coord, sphere_info(j,step_t)%coord, boxsize)< r_LU2 &
  !       ) THEN 
  !       wat_neighbors_number(i) = wat_neighbors_number(i) + 1
  !       !If the number of neighbors are too large, we neglect the extral neighbors (Here, they are water molecules.)
  !       IF (wat_neighbors_number(i) > num_coord_mol_max) EXIT
  !       wat_neighbors_indices(wat_neighbors_number(i),i) = j ! I have to separate the O, and H indies.
  !       !WRITE(*,*) "ATOM Name (k):", sphere_info(k,step_t)%name
  !       !WRITE(*,*) "ATOM Name (j):", sphere_info(j,step_t)%name
  !       !WRITE(*,*) "Atom id (j):", sphere_info(j,step_t)%id
  !       !WRITE(*,*) "distance (k,j):", SQRT(distance2(sphere_info(k,step_t)%coord, sphere_info(j,step_t)%coord, boxsize))
  !    END IF
  !    END DO Inner
  !    !WRITE(*,*) "wat_neighbors_number = ", wat_neighbors_number(i)
  !    !WRITE(*,*) "wat_neighbors_indices:", wat_neighbors_indices(:,i)
  !  END DO Outer
  !END SUBROUTINE make_wat_neighbors_for_X_shell

  CHARACTER(LEN=200) FUNCTION name_list_OH_at_shell_interface(filename,core_atom_str)
    IMPLICIT NONE

    CHARACTER(len=*),INTENT(IN) :: filename 
    character(LEN=*),INTENT(IN) :: core_atom_str        ! For example, core_atom_str can be "I","Na", etc.

    name_list_OH_at_shell_interface=trim(filename)//"_"//core_atom_str//'_OH_at_shell_interface_list.dat'
  END FUNCTION

  CHARACTER(LEN=200) FUNCTION name_c2_at_shell_interface(filename,core_atom_str)
    IMPLICIT NONE

    CHARACTER(len=*),INTENT(IN) :: filename 
    character(LEN=*),INTENT(IN) :: core_atom_str        ! For example, core_atom_str can be "I","Na", etc.

    name_c2_at_shell_interface=trim(filename)//"_"//core_atom_str//'_c2_at_shell_interface.dat'
  END FUNCTION
