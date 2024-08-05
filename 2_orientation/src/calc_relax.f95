!!=======================================================================
!  2021-3-21
!  2024-8-5
!  Author: Huang Gang
!  subroutine calc_relax.f95
!  Autocorrelation of relative positon vector:
!  Auto correlation fuction of auto correlation fucntion of O-H vector.
!  In this function, we calculate the relaxiation of O-H bonds. 
!=======================================================================

  SUBROUTINE calc_relax(boxsize,filename,pos_filename,list_filename,delta_t,natoms,nmovie,&
            n_grid,divx,divy,divz,nb_divx,nb_divy,nb_divz,thickness,surf_info_fortran) 
    !=========================================================
    !Purpose: to obtain C2(t)= P2(u(0)u(t)) for interfacial HB
    !1)the correlation function is the average over N pairs of molecules. 
    !2)One can choose any time step to start, and end at any step for 
    !calculation, instead of the original time step.      
    !3)considered the PBC in an easier way.
    !
    !USE file_module, ONLY: calc_num_of_lines_of_file
    use module_ihb, ONLY: get_total_number_of_lines, &
                     hydrogen_ndx_list_for_oxygen, &
                     distance2, diff_axis, &
                     pm_adh, pm_ahd, & 
                     grid_index, oh_in_surf1, oh_in_surf2,&
                     str, nth, atom, P2
    implicit none
  
    !Declaration
    integer,parameter :: rk=8
    integer, parameter :: d_len=3 ! for storing the length of the character which represents the thickness of the interface

    real(kind=rk),parameter :: rate=0.8d0       ! Condition for cutting off autocorrelation functions
    REAL(KIND=rk),PARAMETER :: max_time_for_corr = 12.0 ! Unit: ps.

    character(LEN=*),INTENT(IN) :: list_filename
    character(LEN=LEN(list_filename)) :: list_filename_temp
    character(LEN=200),INTENT(IN) :: filename, pos_filename
    real(kind=rk),INTENT(IN) :: delta_t      ! ps  Use the time step (delta_t) as a general variable, instead of a parameter.  
    INTEGER,INTENT(IN) :: natoms
    INTEGER,INTENT(IN) :: nmovie
    INTEGER :: nmo_effective, start_step, num_start_points
    real(kind=rk), dimension(3), INTENT(IN) :: boxsize
    integer :: i,iatom,imovie,j,mt,nts
    integer :: nbonds ! number of OH bonds 
    integer, dimension(:,:), allocatable :: idx_bonds
    double precision,dimension(:,:), allocatable :: x,y,z,rx,ry,rz
    double precision,dimension(:), allocatable :: temp_leng ! length of each OH bond at all times
    double precision, dimension(:),allocatable :: corr
    character(len=3), dimension(:), allocatable :: atom_type
    INTEGER, INTENT(IN) :: nb_divx, nb_divy, nb_divz, n_grid 
    REAL(kind=rk), INTENT(IN) :: divx, divy, divz
    REAL(KIND=rk),DIMENSION(2,nmovie,n_grid),INTENT(IN) :: &
        surf_info_fortran

    REAL(kind=rk), INTENT(IN) :: thickness ! the thickness of the instantaneous interfaces
    CHARACTER(len=d_len) :: char_thickness ! for saving the thickness in the files' names
    INTEGER :: ierror 
    integer, allocatable,dimension(:) :: ndx_1
    !==============
    !Initialization
    !==============
    i=0; j=0
    iatom=0
    imovie=0
    mt=0
    nts=0 
    nbonds=0
    nmo_effective = 0
    start_step = 1
    num_start_points = 0

    nmo_effective = nint(max_time_for_corr/delta_t) + 1
    start_step = nint((nmo_effective-1)/5.0) ! Start step of sliding window. Over-using rate is 1 - 1/5 =     4/5
    num_start_points = (nmovie - nmo_effective-1)/start_step + 1
    !Write delta_t
    !write(*,*) "RELAX delta_t = ", delta_t , "(ps)" 
    ! I want to calculate 'nbonds' from the list_filename! Then I need to create an explicit interface, therefore, I defined a function calc_num_of_lines_of_file and put it in a module.
    ! Get the total number of OH bond pairs, i.e., the total number of lines in the list file
    list_filename_temp = list_filename
    nbonds=get_total_number_of_lines(list_filename_temp)

    allocate(ndx_1(nbonds))          
    !============================
    !read data from the list file
    !============================
    !Allocation
    allocate (x(natoms,nmovie),y(natoms,nmovie),z(natoms,nmovie),&
        atom_type(natoms),corr(nmovie) )
    allocate (idx_bonds(2,nbonds),&
        rx(nbonds,nmovie),ry(nbonds,nmovie),rz(nbonds,nmovie) )
    allocate (temp_leng(nmovie))
    !Initialize corr array
    corr = 0.d0

    open(10,file=list_filename_temp,iostat=ierror)
      !Filling the idx_bonds array with the idx of the atoms
      !making the bonds 
      do i=1,nbonds
        read(10,*)idx_bonds(1,i),idx_bonds(2,i)
      end do
    close(10)

    ! read in TRAJECTORY/VELOCITY file in xyz format 
    open(10,file=pos_filename,iostat=ierror)

    do imovie = 1,nmovie
       read(10,*)
       read(10,*)
       do iatom = 1,natoms
          read(10,*)atom_type(iatom),&
               x(iatom,imovie),y(iatom,imovie),z(iatom,imovie)
       enddo
    enddo
    close(10)
    write(6,*)'end of TRAJECTORY reading'

    !First: Calculation for ALL the times of ALL the studied bonds
    !TODO: remove the timeloop.
    do i = 1,nbonds 
       timeloop:DO j = 1, nmovie 
       rx(i,j) = diff_axis(x(idx_bonds(2,i),j),x(idx_bonds(1,i),j),boxsize(1))
       ry(i,j) = diff_axis(y(idx_bonds(2,i),j),y(idx_bonds(1,i),j),boxsize(2))
       rz(i,j) = diff_axis(z(idx_bonds(2,i),j),z(idx_bonds(1,i),j),boxsize(3))
       temp_leng(j) = SQRT(rx(i,j)**2 + ry(i,j)**2 + rz(i,j)**2)
       ! transfer to unit vector 
       rx(i,j) = rx(i,j)/temp_leng(j)
       ry(i,j) = ry(i,j)/temp_leng(j)
       rz(i,j) = rz(i,j)/temp_leng(j)
       END DO timeloop
    enddo
    
    do mt=0, nmo_effective - 1 ! Loop of Dt (=t-t0)
       do nts = 1, nmovie - nmo_effective, start_step ! Loop of t0
          do i = 1, nbonds ! O(5) in wat1
             corr(mt+1) = corr(mt+1) +          &
                  P2(( rx(i,nts)*rx(i,nts+mt) + &
                  ry(i,nts)*ry(i,nts+mt) +      &
                  rz(i,nts)*rz(i,nts+mt)        &
                  )) !Scalar product of vectors
          end do
       end do
       corr(mt+1) = corr(mt+1)/dble(nbonds * num_start_points) ! Normalization by the number of time steps
    end do
    !Normalization by corr(1)=<r(0).r(0)>
    do mt=nmo_effective, 1, -1
      corr(mt)=corr(mt)/corr(1)
    enddo

    char_thickness = nth(str(thickness),d_len)
    open(10,file=trim(filename)//'_c2_ihb_' &
      //char_thickness//'.dat')
    do i = 1, nmo_effective
       write(10,*)(i-1)*delta_t,corr(i)
    enddo
    write(6,*)'correlation results written in file '//trim(filename)//'_c2_ihb_'//char_thickness//'.dat'
    close(10)

    !============
    !deallocation
    !============
    deallocate (x,y,z,rx,ry,rz,atom_type,idx_bonds,corr,temp_leng)

  END SUBROUTINE
