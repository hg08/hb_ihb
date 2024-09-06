      SUBROUTINE gfreeoh(filename,pos_filename,list,list_H,natoms)
      !
      !Purpose: To generate the index list of all Hydrogen atoms 
      !============
      ! 2024/8
      ! huang gang
      !============
      USE parameter_shared, ONLY: ndx_O, ndx_H
      use water_molecule_types

      implicit none
      !========================
      !Parameters and variables
      !========================
      character(LEN=200),INTENT(IN) :: filename            ! specific filename to analyz data
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      character(LEN=200),INTENT(INOUT) :: list
      character(LEN=200),INTENT(INOUT) :: list_H
      INTEGER, INTENT(IN) :: natoms

      ! Local variables
      integer,parameter :: rk=8              
      real,parameter :: r_ohc=1.44             ! rOH (1.2**2)
      integer :: i,j,nmovie,iatom,& 
                    m1,m2,&
                    i1,i2,ii,jj,i_O,i_H
      real,allocatable,dimension (:,:) :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      !integer,allocatable,dimension (:) :: ndx_O,ndx_H

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
      read(10,*) ! Neglect the first two lines
      read(10,*)                  
      i=0
      ii=0
      do iatom= 1,natoms
          read (10,*)atom_type(iatom),x(iatom,nmovie),& 
                   y(iatom,nmovie),z(iatom,nmovie)
          if (trim(atom_type(iatom)) .eq. 'O') then
              i=i+1
          elseif(trim(atom_type(iatom)) .eq. 'H') then
              ii=ii+1
          endif
      enddo
      i_O=i
      i_H=ii
      close(10)
      !=======================================================
      ! Calculate the number of O (H) atoms in water molecules
      !=======================================================
      allocate(ndx_O(i_O,2)) ! this should be put after i_O is defined
      allocate(ndx_H(i_H,2)) ! this should be put after i_H is defined
      i=0
      ii=0
      jj=0
      do iatom=1,natoms
          if (trim(atom_type(iatom)) .eq. 'O')then
                 i=i+1      
                 ndx_O(i,1)=i ! O index in all O atoms.
                 ndx_O(i,2)=iatom ! O index in all atoms.
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ii=ii+1
                 ndx_H(ii,1)=ii ! H index in all H atoms.
                 ndx_H(ii,2)=iatom ! H index in all atoms. 
          else
          endif
      enddo 
      deallocate(atom_type)
      !========================
      !Producing the O-O-H list      
      !========================
      list=trim(filename)//'_O_list.dat'
      open(20,file=list)
          do i1=1, i_O-1 ! No O atom can not be bonded to itself 
              !m1=ndx_O(i1,2)
              do i2=i1+1, i_O ! Acceptor 
                  !m2=ndx_O(i2,2)
                  write(20,*)  ndx_O(i1,1), ndx_O(i1,2), ndx_O(i2,1), ndx_O(i2,2)
                  !write(20,*) m1, m2
              enddo
          enddo
      close(20)
      !====================
      !Producing the H list      
      !====================
      list=trim(filename)//'_H_list.dat'
      open(20,file=list)
          do i1=1, i_H 
              write(20,*)  ndx_H(i1,1), ndx_H(i1,2)
          enddo
      close(20)
      !deallocate(ndx_O,ndx_H,x,y,z)
      deallocate(x,y,z)
      !====================
      END SUBROUTINE
