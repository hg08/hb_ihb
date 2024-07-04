      SUBROUTINE ghbacf_interface_n_pbc_format2(boxsize, delta_t0, &
          filename, pos_filename,list_filename,n_samples,nat, ns, &
          criterion)
      !2020-5-29: simplifying the definition of hb (without r13)
      !2020-9-18: In this function, we define interfacial HB. 
      !Calculate the correlation n(t) of interfacial HBs.
      !====================================================================
      !Purpose: to obtain n_HB(t) for interfacial HBs.
      !1)the correlation function is the average over N pairs of molecules. 
      !2)One can choose any time step to start, and end at any step for 
      !calculation, instead of the original time step.      
      !3)considered the PBC in an easier way.
      !4) surf_filename have to be provided, which include the surface's 
      ! coordinates. In the current version, the surf trajectory is 
      ! obtained
      !5)A local parameter thickness is needed to define the thickness of 
      ! the instantaneous interfaces for the system.
      !====================================================================
      ! INPUT: 
      ! box size
      ! time step
      ! name of system
      ! name of trajectory
      ! name of list
      ! nmo
      ! nat
      ! Hydrogen bond definition (1 or 2)
      ! surf_filename
      !=========================================================
      !Modules
      use tools, ONLY: get_nwat => get_total_numer_of_lines, &
                       h_ndx_list => hydrogen_ndx_list, &
                       dist2 => distance2, &
                       pmADH, pmAHD, &
                       grid_index, pair_in_surf1, pair_in_surf2,&
                       str, nth
      USE atom_module
      USE parameter_shared
      USE traj_format2
      USE surf_traj
      USE surf_module, ONLY: surf_info

      implicit none
      !==========
      !parameters
      !==========
      character(LEN=200),INTENT(INOUT) :: filename,pos_filename
      character(LEN=200),INTENT(IN) :: list_filename
      integer,INTENT(IN) :: criterion
      INTEGER, INTENT(IN) :: nat ! number of atoms
      INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
      integer,parameter :: rk=8 ! local 
      real(kind=rk), dimension(3), INTENT(IN) :: boxsize
    
      !Local variables
      real,parameter :: rooc=12.25                 ! cutoff distance of rOO (3.5**2 )
      real,parameter :: cosPhiC123=0.866           ! 1.732/2; phiC=pi/6.
      real,parameter :: cosPhiC132=-0.5            ! -1./2; phiC132=2pi/3.
      real(kind=rk),parameter :: h_min=0.5  ! condition for the existence of a h-bond for a step
      real(kind=rk),parameter :: hb_min=0.5 ! condition for the existence of h-bond for a pair of water molecules
      real(kind=rk) :: r13, cosphi, pm, cosphi_, pm_
      real(kind=rk) :: r21, r31, r32, r23 ! For the second criterion of HB
      real(kind=rk) ::qj, tot_hb, delta_t, delta_t0, hb_per_frame, ave_h
      real(kind=rk), dimension(3) :: r1, r2, r3 ! pbc 
      integer :: m1,m2,m3,mt,nqj,tot_nhb,n_bonded_pairs,ns
      real(kind=rk),allocatable,dimension (:)  :: h,h_d,hb,corr_n
      real,allocatable,dimension (:,:)         :: x,y,z
      integer,allocatable,dimension(:)         :: ndx_1,ndx_2,nhb_exist
      integer,dimension(4)   :: ndx_3_list
      real(kind=rk)  :: scalar 
      logical,allocatable,dimension (:)  :: hb_exist
      INTEGER  :: nmo  ! nmo is not necessary, we set nmo = n_samples, because we do not want to change too much
      INTEGER :: nwat ! number of water molecules
      INTEGER :: i,j,k,jj 
      CHARACTER(len=1) :: char_thickness ! for saving the thickness in the files' names
      INTEGER :: index_mol1, index_mol2
      LOGICAL :: condition1,condition2
      !==============
      !Initialization
      !==============
      ave_h =0.0; scalar =0.0
      pm =0.0; cosphi =0.0
      r21 = 0.0; r23 = 0.0
      r31 = 0.0; r13= 0.0; r32 = 0.0
      hb_per_frame = 0.0; tot_hb = 0.0
      r1 = 0.0; r2 = 0.0; r3 = 0.0
      nmo = n_samples; nwat=0 
      ndx_3_list=0
      index_mol1=0; index_mol2=0
      condition1=.FALSE.
      condition2=.FALSE.

      !To obtain the total number of water pairs
      nwat=get_nwat(list_filename)
      write(*,*) 'ghbacf_c: # of water pairs (nwat) =', nwat
      allocate(ndx_1(nwat))          
      allocate(ndx_2(nwat))          
      !============================
      !read data from the list file
      !============================
      open(10,file=list_filename)     
      do k=1,nwat
          read(10,*)ndx_1(k),ndx_2(k)
      enddo
      close(10)
      !============================

      delta_t=ns*delta_t0  ! unit: ps
      write(*,*) "New total steps (nmo):", nmo
      allocate(x(nat,nmo))
      allocate(y(nat,nmo))
      allocate(z(nat,nmo))
      allocate(h(nmo))
      allocate(h_d(nmo))
      allocate(hb(nwat)) ! Average H-bonded population 
      allocate(nhb_exist(nwat))
      !============================
      !read in surf trajectory file 
      !============================
      open(20,file='surf_traj.dat')
           
      !==============================      
      ! Calculate correlation n_HB(t)  
      !==============================      
      allocate(corr_n(nmo))
      allocate(hb_exist(nmo))
      corr_n(:)=0.0      
      tot_hb=0.0
      tot_nhb=0
      h=0.0; h_d=0.0
      hb(:)=0.0
      nhb_exist(:)=0 
      !=============
      !The main loop
      !=============      
      kLOOP: do k=1,nwat
        qj=0
        nqj=0 ! The number of bonded times for k-th form of quasi-HB 
        m1=ndx_1(k)
        m2=ndx_2(k)
        ndx_3_list= & 
            h_ndx_list(ndx_1(k),ndx_2(k),pos_filename,nat,boxsize)
        write(*,*) "The ",k,"-th pair: ndx_of H (1st,2nd,3rd,4th):",& 
            ndx_3_list(1), ndx_3_list(2), ndx_3_list(3), ndx_3_list(4)
        ! Calculate h_d(j)
        TIME: do jj =1, nmo
          h(jj)=0.0
          hb_exist(jj)=.False.

          ! Check if the pairs are located in one of the interfaces 
          index_mol1 =grid_index(atom_info(m1,jj)%coord(1), &
              atom_info(m1,jj)%coord(2),divx,divy,nb_divx) 
          WRITE(*,*) "index_mol1 = ",index_mol1
          index_mol2 =grid_index(atom_info(m2,jj)%coord(1), &
              atom_info(m2,jj)%coord(2),divx,divy,nb_divx) 

          !For surf 1
          condition1 = pair_in_surf1(surf_info(index_mol1,jj)%coord(1),&
              atom_info(m1,jj)%coord(3), &
              surf_info(index_mol2,jj)%coord(1), &
              atom_info(m2,jj)%coord(3),thickness ) 

          !For surf 2 
          condition2 = pair_in_surf2(surf_info(index_mol1,jj)%coord(2),&
              atom_info(m1,jj)%coord(3), &
              surf_info(index_mol2,jj)%coord(2), &
              atom_info(m2,jj)%coord(3),thickness ) 
          WRITE(*,*) condition1, condition2 

          !This condition is the additional condition for the establishment 
          ! of interface hydrogen bonds, which is the core of this method. 
          IF (condition1 .OR. condition2) THEN

              ! A LOOP on ndx_3_list
              HYDROGEN: DO j=1,4
                  !Try each Hydrogen bond, there are totally 4 hydrogen atoms, for a water-water pair
                  m3=ndx_3_list(j)
                  r1=(/atom_info(m1,jj)%coord(1),atom_info(m1,jj)%coord(2),&
                         atom_info(m1,jj)%coord(3) /)
                  r2=(/atom_info(m2,jj)%coord(1),atom_info(m2,jj)%coord(2),&
                         atom_info(m2,jj)%coord(3) /)
                  r3=(/atom_info(m3,jj)%coord(1),atom_info(m3,jj)%coord(2),&
                         atom_info(m3,jj)%coord(3) /)
                  r21 = dist2(r1, r2, boxsize)

                  if (criterion==1) then
                      !Follow the first cirterion of HB.
                      r23 = dist2(r3, r2, boxsize) 
                      r13 = dist2(r3, r1, boxsize) 
                      pm = pmADH(r1,r2,r3,boxsize)  ! if H is bound to O2
                      pm_ = pmADH(r2,r1,r3,boxsize) ! if H is bound to O1
                      cosphi= pm/(sqrt(r21*r23))
                      cosphi_= pm_/(sqrt(r21*r23))
                      IF (r21 .lt. rooc ) then
                          h_d(jj)=1.0 
                      ENDIF
                      IF ( (r21 .lt. rooc).and.((cosphi.gt.cosPhiC123) .or. &
                          (cosphi_ .gt. cosPhiC123) )        &
                         ) THEN
                          h(jj)=1.0 
                          hb_exist(jj) = .True.
                          qj=qj+h(jj) ! To calculate ave population of HB over all starting points for one pair of water                           
                          nqj=nqj+1
                          EXIT ! if we know that two pair of molecule is bonded at frame jj, then we go to check the next frame (jj+1)
                      ENDIF 
                  elseif (criterion==2) then
                      !Follow the second cirterion of HB.
                      r31 = dist2(r1,r3,boxsize) 
                      r32 = dist2(r2,r3, boxsize) 
                      pm= pmAHD(r1,r2,r3,boxsize)
                      cosphi= pm/(sqrt(r31*r32))
                      IF (r21 .lt. rooc ) then
                          h_d(jj)=1.0 
                      ENDIF
                      IF ((r21 .lt. rooc ).and.(cosphi .lt. cosPhiC132) &
                         ) THEN
                          h(jj)=1.0 
                          hb_exist(jj) = .True.
                          qj=qj+h(jj) ! To calculate ave population of HB over all starting points for one pair of water                           
                          nqj=nqj+1
                          EXIT ! if we know that two pair of molecule is bonded at frame jj, then we exit and go to check the next frame (jj+1)
                      ENDIF 
                   endif
               END DO HYDROGEN
           ENDIF
        enddo TIME 
        hb(k)=qj 
        nhb_exist(k)=nqj
        tot_hb=tot_hb + qj
        tot_nhb=tot_nhb+nhb_exist(k)
        !==========================================
        !Calcualte the correlation function n_HB(t)
        !==========================================
        if (hb(k)>hb_min) then
            do mt=0,nmo-1    ! time interval
                scalar=0.d0
                !do j=1,nmo-mt-1
                do j=1,nmo-mt
                    scalar=scalar+h(j)*(1-h(j+mt))*h_d(j+mt)  
                enddo
                corr_n(mt+1)=corr_n(mt+1)+scalar !sum_C_k(t)
            enddo
        endif
      enddo kLOOP   
      deallocate(hb_exist,nhb_exist,h_d)
      hb_per_frame = tot_hb/nmo
      write(*,*) "Total H-bonds exists in history: ", tot_hb
      ave_h = hb_per_frame/nwat 
      write(*,*) "hb per frame:", hb_per_frame
      write(*,*) "<h> =", ave_h
      !=========================================
      !Calculate the number of ever bonded pairs
      !=========================================
      n_bonded_pairs=0 
      do k=1,nwat
          if (hb(k)>hb_min) then
              n_bonded_pairs=n_bonded_pairs+1      
          endif
      enddo
      !==============================
      !Normalization of n_HB(t) step1
      !==============================
      do mt=0,nmo-1! time interval
          corr_n(mt+1)=corr_n(mt+1)/(nmo-mt)
      enddo
      corr_n = corr_n /nwat
      !Normalization step2
      corr_n = corr_n/ave_h
      deallocate(x,y,z,ndx_1,ndx_2)          
     !====================================
     !Write the correlation
     !n_HB(t) for the interfacial HB (ihb)     
     !====================================
      char_thickness = nth(str(nint(thickness)),1)
      open(10,file=trim(filename)//'_wat_pair_hbacf_n_ihb_'//&
        char_thickness//'.dat')
        do i=1,nmo
            write(10,*)(i-1)*delta_t,corr_n(i)
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_wat_pair_hbacf_n_ihb_'//char_thickness//'.dat'
      close(10)
     !=====================
     !Write the correlation
     !ln(n_HB(t))     
     !=====================
      open(10,file=trim(filename)//'_wat_pair_hbacf_log_n_ihb_'//&
        char_thickness//'.dat')
        do i=1,nmo
            write(10,*)(i-1)*delta_t,log(corr_n(i))
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_wat_pair_hbacf_log_n_ibh_'//char_thickness//'.dat'
      close(10)
     !===========
     ! Print <h>      
     !===========      
      open(10,file=trim(filename)//'_wat_pair_ave_h_by_n_ihb_'//&
        char_thickness//'.dat')
        write(10,*) 'Ave. No. bonds:', hb_per_frame
        write(10,*) '<h>:',ave_h
        write(6,*)'written in '//trim(filename)//&
                  '_wat_pair_ave_h_by_n_ihb_'//char_thickness//'.dat'
      close(10)
      deallocate (h,corr_n,hb)
      END SUBROUTINE 
