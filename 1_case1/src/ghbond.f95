      SUBROUTINE ghbond(filename,pos_filename,list,natoms)
          !
          !Purpose: To generate the index list of all Oxygen-Oxygen pairs (quasi-bonds)
          !
          !==================
          !Date: 2020/07
          !Author: Huang Gang
          !==================
          implicit none
          !========================
          !Parameters and variables
          !========================
          character(LEN=200),INTENT(IN) :: filename ! specific filename to analyz data
          character(LEN=200),INTENT(INOUT) :: list
          INTEGER, INTENT(IN) :: natoms
          character(LEN=200),INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data

          ! Local variables
          integer :: i,ii,ierror,i1,i2,i_O,i_H,iatom,imovie,&
                     j,jj,m1,m2,nmovie 
          character(LEN=3),allocatable,dimension (:) :: atom_type
          integer,allocatable,dimension (:) :: ndx_H,ndx_O
          real,allocatable,dimension (:,:) :: x,y,z

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
          open(10,file=trim(pos_filename),STATUS='OLD',ACTION='READ',IOSTAT=ierror)  
          do imovie=1,nmovie
              read(10,*)!Neglect data of this line
              read(10,*)                  
              i=0
              ii=0
              do iatom= 1,natoms
                  read (10,*)atom_type(iatom),x(iatom,imovie),& 
                           y(iatom,imovie),z(iatom,imovie)
                  if (trim(atom_type(iatom)) .eq. 'O') then
                      i=i+1
                  elseif(trim(atom_type(iatom)) .eq. 'H') then
                      ii=ii+1
                  endif
              enddo
          enddo
          i_O=i
          i_H=ii
          write(6,*)'Number of O','    Number of H'
          write(6,*)i_O,i_H
          close(10)
          !=======================================================
          ! Calculate the number of O (H) atoms in water molecules
          !=======================================================
          allocate(ndx_O(i_O)) ! this should be put after i_O is defined
          allocate(ndx_H(i_H)) ! this should be put after i_H is defined
          i=0; ii=0; jj=0
          do iatom=1,natoms
              if (trim(atom_type(iatom)) .eq. 'O')then
                     i=i+1      
                     ndx_O(i)=iatom
              elseif(trim(atom_type(iatom)) .eq. 'H') then
                     ii=ii+1
                     ndx_H(ii)=iatom
              else
              endif
          enddo 
          deallocate(atom_type)
          !========================
          !Producing the O-O-H list      
          !========================
          list=trim(filename)//'_O_list.dat'
          open(20,file=list,STATUS='REPLACE',ACTION='WRITE')
              do i1=1, i_O-1! No O atom can not be bonded to itself 
                  m1=ndx_O(i1)
                  do i2=i1+1, i_O 
                      m2=ndx_O(i2)
                      do j =1, 1
                          write(20,*) m1,m2
                      enddo
                  enddo
              enddo
          close(20)
          deallocate(ndx_O,ndx_H,x,y,z)
      
      END SUBROUTINE ghbond 
