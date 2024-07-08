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
          CHARACTER(LEN=200),INTENT(IN) :: filename ! specific filename to analyz data
          CHARACTER(LEN=200),INTENT(INOUT) :: list
          INTEGER, INTENT(IN) :: natoms
          CHARACTER(LEN=200),INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data

          ! Local variables
          INTEGER :: i,ii,ierror,i1,i2,i_O,i_H,iatom,imovie,&
                     j,jj,m1,m2,nmovie 
          CHARACTER(LEN=3),ALLOCATABLE,dimension (:) :: atom_type
          INTEGER,ALLOCATABLE,dimension (:) :: ndx_H,ndx_O
          REAL,ALLOCATABLE,dimension (:,:) :: x,y,z

          !==================
          !Read data in input
          !==================
          nmovie=1!Number of movie steps reqired.
          ALLOCATE(atom_type(natoms))
          ALLOCATE(x(natoms,nmovie))
          ALLOCATE(y(natoms,nmovie))
          ALLOCATE(z(natoms,nmovie))
          !=======================
          !read in trajectory file 
          !=======================
          OPEN(10,file=trim(pos_filename),STATUS='OLD',ACTION='READ',IOSTAT=ierror)  
          DO imovie=1,nmovie
              READ(10,*)!Neglect data of this line
              READ(10,*)                  
              i=0
              ii=0
              DO iatom= 1,natoms
                  READ (10,*)atom_type(iatom),x(iatom,imovie),& 
                           y(iatom,imovie),z(iatom,imovie)
                  IF (trim(atom_type(iatom)) .eq. 'O') THEN
                      i=i+1
                  ELSEIF(trim(atom_type(iatom)) .eq. 'H') THEN
                      ii=ii+1
                  ENDIF
              ENDDO
          ENDDO
          i_O=i
          i_H=ii
          WRITE(6,*)'Number of O','    Number of H'
          WRITE(6,*)i_O,i_H
          CLOSE(10)
          !=======================================================
          ! Calculate the number of O (H) atoms in water molecules
          !=======================================================
          ALLOCATE(ndx_O(i_O)) ! this should be put after i_O is defined
          ALLOCATE(ndx_H(i_H)) ! this should be put after i_H is defined
          i=0; ii=0; jj=0
          DO iatom=1,natoms
              IF (trim(atom_type(iatom)) .eq. 'O')THEN
                     i=i+1      
                     ndx_O(i)=iatom
              ELSEIF(trim(atom_type(iatom)) .eq. 'H') THEN
                     ii=ii+1
                     ndx_H(ii)=iatom
              ELSE
              ENDIF
          ENDDO 
          deALLOCATE(atom_type)
          !========================
          !Producing the O-O-H list      
          !========================
          list=trim(filename)//'_O_list.dat'
          OPEN(20,file=list,STATUS='REPLACE',ACTION='WRITE')
              DO i1=1, i_O-1! No O atom can not be bonded to itself 
                  m1=ndx_O(i1)
                  DO i2=i1+1, i_O 
                      m2=ndx_O(i2)
                      DO j =1, 1
                          WRITE(20,*) m1,m2
                      ENDDO
                  ENDDO
              ENDDO
          CLOSE(20)
          deALLOCATE(ndx_O,ndx_H,x,y,z)
      
      END SUBROUTINE ghbond 
