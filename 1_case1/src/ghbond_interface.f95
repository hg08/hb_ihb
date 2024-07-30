      SUBROUTINE ghbond_interface(filename,list,n_samples,nat,arr1,arr2)
      !
      !Purpose: To generate the index list of all Oxygen-Oxygen pairs (quasi-bonds) which are once located in
      !         the neat water interface. There's an inportant idea in this code: We assume that ergodity is not satisfied,
      !         therefore, for a given layer thickness, we list all O-O pairs' indices in one list file, 
      !         and calculate the hydrogen bond dynamics for all these pairs (quasi-HBs), the correlation of some
      !         pairs may calculated serveral times. It is exactly what we want, because the more the pair is
      !         considered, the more weights this pair will have in the statistis. After all the correlations
      !         are calculated, we calculate the average correlation over all these correlations.
      !Output: A list file of O-O pairs. Note that one O-O pair may occurs more than once, because the O-O
      !        pair may in interface at both time t_i and t_j.
      !============
      ! 2020/10
      ! huang gang
      !============

      IMPLICIT NONE
      !========================
      !Parameters and variables
      !========================
      CHARACTER(LEN=200),INTENT(IN) :: filename ! specific filename to analyz data
      CHARACTER(LEN=200),INTENT(INOUT) :: list
      INTEGER,INTENT(IN) :: n_samples ! number of samples
      !INTEGER,INTENT(IN) :: nmo ! steps of trajectory
      INTEGER,INTENT(IN) :: nat ! number of atoms
      !To save the indices of the molecules for generating list file, we define an array for each time point (jj, in this code)
      INTEGER,DIMENSION(n_samples,nat),INTENT(INOUT) :: arr1, arr2 ! Jie: updated to array 1 and 2.

      ! Local variables
      INTEGER,PARAMETER :: step=100 ! The parameter 'step' should not be too small, otherwise, you will waste your time on many repeated calculation. 
                                    ! Here, 100 means : "We select the molecules in interface every 100*ns steps."
      INTEGER :: i,n,m1,m2,i1,i2,jj
      INTEGER,DIMENSION (nat) :: ndx_O
     
      !Initialization
      ndx_O=0; n=0

      !====================================================
      !Producing the O-O list for O atom pairs in interface     
      !====================================================
      list=trim(filename)//'_O_in_interface_list.dat'
      ! DO NOT USE "status='NEW' ".
      OPEN(20,file=list, STATUS='REPLACE', ACTION='WRITE') ! Set the status to 'REPLACE', otherwise, the compiler can not find this file.
        write(*,*) 'Debug: open list file'
        !DO-LOOP on each row of the indx array 'arr1'
        ROW: DO jj=1,n_samples,step  ! start, end [, step]
          ndx_O(:)=arr1(jj,:) 
          !===============================================
          !Calculate the number of O atoms in this time jj
          !===============================================
          n=0
          DO i=1,nat
            IF (ndx_O(i)>0) THEN
               n=n+1
            ENDIF
          ENDDO
          DO i1=1,n-1 ! No O atom can not be bonded to itself 
            m1=ndx_O(i1)
            DO i2=i1+1,n 
              m2=ndx_O(i2)
              WRITE(20,*) m1,m2
            ENDDO
          ENDDO
        ENDDO ROW
        
        !DO-LOOP on each row of the indx array 'arr2'
        ROW2: DO jj=1,n_samples,step  ! start, end [, step]
          ndx_O(:)=arr2(jj,:) 
          !===============================================
          !Calculate the number of O atoms in this time jj
          !===============================================
          n=0
          DO i=1,nat
            IF (ndx_O(i)>0) THEN
               n=n+1
            ENDIF
          ENDDO

          DO i1=1,n-1 
            m1=ndx_O(i1)
            DO i2=i1+1,n 
              m2=ndx_O(i2)
              WRITE(20,*) m1,m2
            ENDDO
          ENDDO
        ENDDO ROW2

      CLOSE(20)
      !====================
      END SUBROUTINE ghbond_interface 
