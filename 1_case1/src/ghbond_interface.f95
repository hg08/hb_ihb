      SUBROUTINE ghbond_interface(filename,list,nmo,nat,arr)
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

      implicit none
      !========================
      !Parameters and variables
      !========================
      character(LEN=200),INTENT(IN) :: filename            ! specific filename to analyz data
      character(LEN=200),INTENT(INOUT) :: list
      INTEGER,INTENT(IN) :: nmo ! steps of trajectory
      INTEGER,INTENT(IN) :: nat ! number of atoms
      !To save the indices of the molecules for generating list file, we define an array for each time point (jj, in this code)
      INTEGER,DIMENSION(nmo,nat),INTENT(INOUT) :: arr

      ! Local variables
      INTEGER,PARAMETER :: step=100 ! The parameter 'step' should not be too small, otherwise, you will waste your time on many repeated calculation. 
                                    ! Here, 100 means : "We select the molecules in interface every 100*ns steps."
      integer :: i,n,m1,m2,i1,i2,jj
      integer,dimension (nat) :: ndx_O
     
      !Initialization
      ndx_O=0;n=0

      !====================================================
      !Producing the O-O list for O atom pairs in interface     
      !====================================================
      list=trim(filename)//'_O_in_interface_list.dat'
      ! DO NOT USE "status='NEW' ".
      open(20,file=list, STATUS='REPLACE', ACTION='WRITE') ! Set the status to 'REPLACE', otherwise, the compiler can not find this file.
      
      !DO-LOOP on each row of the indx array 'arr'
      ROW: DO jj=1,nmo,step  ! start, end [, step]
        ndx_O(:)=arr(jj,:) 
        !===============================================
        !Calculate the number of O atoms in this time jj
        !===============================================
        n=0
        do i=1,nat
          if (ndx_O(i)>0) then
             n=n+1
          endif
        enddo

        do i1=1,n-1 ! No O atom can not be bonded to itself 
          m1=ndx_O(i1)
          do i2=i1+1,n 
            m2=ndx_O(i2)
            write(20,*) m1,m2
          enddo
        enddo
      ENDDO ROW

      close(20)
      !====================
      END SUBROUTINE ghbond_interface 
