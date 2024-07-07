      SUBROUTINE oh_interface(boxsize,filename,pos_filename,list,nmo,nat,arr)
      !
      !Purpose: To generate the index list of all Oxygen atoms which are once located in
      !         the neat water interface. 
      !         There's an inportant idea in this code: 
      !         IN each divided piece of the trajectory, we assume that ergodity is not satisfied.
      !         Therefore, for a given layer thickness, we list all O atoms' indices in one list file, 
      !         and calculate C2(t) dynamics for all these OH bonds, the correlation of some
      !         OH bonds may calculated serveral times. It is exactly what we want, because the more the OH's is
      !         considered, the more weights this OH bond will have in the statistis. After all the correlations
      !         are calculated, we calculate the average correlation over all these correlations.
      !Output: A list file of O atoms. Note that one O atom may occurs more than once, because the O atom
      !        may in interface at both time t_i and t_j.
      !============
      ! 2020/10
      ! huang gang
      !============
      USE module_ihb, ONLY: hydrogen_ndx_list_for_oxygen
      implicit none
      !========================
      !Parameters and variables
      !========================
      integer, parameter :: rk=8              
      character(LEN=200),INTENT(IN) :: filename,pos_filename            ! specific filename to analyz data
      character(LEN=200),INTENT(INOUT) :: list
      INTEGER,INTENT(IN) :: nmo ! steps of trajectory
      INTEGER,INTENT(IN) :: nat ! number of atoms
      real(kind=rk), dimension(3), INTENT(IN) :: boxsize
      !To save the indices of the molecules for generating list file, we define an array for each time point jj.
      INTEGER,DIMENSION(nmo,nat),INTENT(INOUT) :: arr

      ! Local variables
      INTEGER,PARAMETER :: step=25 ! The parameter 'step' should not be too small, otherwise, you will waste your time on many repeated calculation. 
                                    ! Here, 100 means : "We select the molecules in interface every 100*ns steps."
      integer :: i,n,jj
      integer,dimension (nat) :: ndx_O
      integer, dimension(2)   :: ndx_3_list
     
      !Initialization
      ndx_O=0;n=0

      !======================================
      !Producing the O atom list in interface     
      !======================================
      list=trim(filename)//'_OH_at_interface_list.dat'
      open(20,file=list)

      !DO-LOOP on each row of the indx array 'arr'
      ROW: DO jj=1,nmo,step  ! start, end [, step]
        ndx_O(:)=arr(jj,:) 
        !==========================================
        !Calculate the number of O atoms in time jj
        !==========================================
        n=0
        do i=1,nat
          if (ndx_O(i)>0) then
            ndx_3_list=hydrogen_ndx_list_for_oxygen(ndx_O(i),pos_filename,nat,boxsize)
            write(20,*) ndx_O(i), ndx_3_list(1)
            write(20,*) ndx_O(i), ndx_3_list(2)
            n=n+2
          endif
        enddo
        write(*,*) "n = ", n
      ENDDO ROW

      close(20)
      !====================
      END SUBROUTINE oh_interface 
