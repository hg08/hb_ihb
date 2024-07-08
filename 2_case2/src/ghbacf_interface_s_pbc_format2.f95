      SUBROUTINE ghbacf_interface_s_pbc_format2(boxsize,delta_t0,&
          filename,pos_filename,list_filename,n_samples,nat,&
          ns,criterion,donor,accepter)

      !2020-5-29: simplifying the definition of hb (without r13)
      !========================================================================
      !Purpose: to obtain S_HB(t): S_{HB}(t)=<h(0)H(t)>/<h> for interfacial HB
      !1)the correlation function is the average over N pairs of molecules. 
      !2)One can choose any time step to start,and end at any step for 
      !calculation,instead of the original time step.      
      !3)considered the PBC in an easier way.
      !========================================================================
      ! INPUT: 
      ! box size
      ! time step
      ! name of system
      ! name of trajectory
      ! name of list
      ! nmo
      ! nat
      ! time step for calculating CORRELATION
      ! Hydrogen bond definition (1 or 2)
      ! surf_filename
      !=========================================================
      !Modules
      use tools,ONLY: get_nwat=> get_total_numer_of_lines,&
                       h_ndx_list=> hydrogen_ndx_list,&
                       dist2=> distance2,&
                       pmADH,pmAHD,& 
                       grid_index,pair_in_surf1,pair_in_surf2,&
                       str,nth
      USE atom_module
      USE parameter_shared
      USE traj_format2
      USE surf_traj
      USE surf_module,ONLY: surf_info
      implicit none
      
      !==========
      !parameters
      !==========
      integer,parameter :: rk=8 ! local 

      character(LEN=*),INTENT(IN) :: accepter 
      character(LEN=*),INTENT(IN) :: donor 
      character(LEN=200),INTENT(INOUT) :: filename,pos_filename
      character(LEN=200),INTENT(IN) :: list_filename
      integer,INTENT(IN) :: criterion
      INTEGER,INTENT(IN) :: nat ! number of atoms
      INTEGER,INTENT(IN) :: n_samples  !n_samples=INT(nmo/ns)
      INTEGER,INTENT(IN) :: ns  
      REAL(kind=rk),dimension(3),INTENT(IN) :: boxsize
      REAL(KIND=rk),INTENT(IN) :: delta_t0 

      !Local variables
      REAL(kind=rk),parameter :: rooc=12.25                 ! cutoff distance of rOO (3.5**2 )
      REAL(kind=rk),parameter :: cosPhiC123=0.866           ! 1.732/2;phiC=pi/6.
      REAL(kind=rk),parameter :: cosPhiC132=-0.5            ! -1./2;phiC132=2pi/3.
      REAL(kind=rk),parameter :: h_min=0.5 ! condition for the existence of a h-bond for a step
      REAL(kind=rk),parameter :: hb_min=0.5 ! condition for the existence of h-bond for a pair of water molecules
      REAL(kind=rk) :: r13,cosphi,pm,cosphi_,pm_,norm_rr
      REAL(kind=rk) :: r21,r31,r32,r23 ! For the second criterion of HB
      REAL(kind=rk) :: qj,tot_hb,delta_t,hb_per_frame,ave_h,hh
      REAL(kind=rk),dimension(3) :: r1,r2,r3 ! pbc 
      integer :: m,m1,m2,m3,mt,nqj,tot_nhb,n_bonded_pairs
      REAL(kind=rk),allocatable,dimension (:)  :: h,hb,corr_hh
      REAL(kind=rk),allocatable,dimension (:,:) :: x,y,z
      integer,allocatable,dimension(:) :: ndx_1,ndx_2,nhb_exist,&
          nhb_start
      integer,dimension(4) :: ndx_3_list
      REAL(kind=rk) :: scalar_hh 
      REAL(kind=rk) :: tau_a ! average lifetime of HBs 
      logical,allocatable,dimension (:)  :: hb_exist,hb_start
      INTEGER :: nmo  ! nmo is not necessary,we set nmo=n_samples,because we do not want to change too much
      INTEGER :: nwat ! number of water molecules
      INTEGER :: i,j,k,jj 
      CHARACTER(len=1) :: char_thickness ! for saving the thickness in the files' names
      CHARACTER(len=1) :: char_criterion ! for saving the type of HB definition: ADH (1) or AHD (2)
      INTEGER :: index_mol1,index_mol2
      LOGICAL :: condition1,condition2

      !==============
      !Initialization
      !==============
      ave_h=0.0;scalar_hh=0.0
      delta_t=0.0
      pm=0.0;cosphi=0.0
      r21=0.0;r23=0.0
      r31=0.0;r13=0.0;r32=0.0
      hb_per_frame=0.0;tot_hb=0.0
      r1=0.0;r2=0.0;r3=0.0
      nmo=n_samples;nwat=0 
      ndx_3_list=0
      norm_rr=0.0
      index_mol1=0;index_mol2=0
      condition1=.FALSE.
      condition2=.FALSE.
      char_criterion=""

      !To obtain the total number of water pairs
      nwat=get_nwat(list_filename)
      WRITE(*,*) 'ghbacf_c: # of water pairs (nwat)=',nwat
      allocate(ndx_1(nwat))          
      allocate(ndx_2(nwat))          
      !============================
      !read data from the list file
      !============================
      OPEN(10,file=list_filename)     
      do k=1,nwat
          read(10,*)ndx_1(k),ndx_2(k)
      enddo
      close(10)
      !============================

      delta_t=ns*delta_t0  ! unit: ps
      WRITE(*,*) "New total steps (nmo):",nmo
      allocate(x(nat,nmo))
      allocate(y(nat,nmo))
      allocate(z(nat,nmo))
      allocate(h(nmo))
      allocate(hb(nwat))    ! Average H-bonded population 
      allocate(nhb_exist(nwat))
      allocate(nhb_start(nwat))
      !====================================
      ! Calculate <h(0)h(t)>/<h>  
      ! Notice here <> is not average over
      ! different pairs of water molecules,
      ! and over all starting time points i
      ! with h(i)=1.
      !====================================     
      allocate(corr_hh(nmo))
      allocate(hb_exist(nmo))
      allocate(hb_start(nmo))
      ! loop
      corr_hh(:)=0.0      
      tot_hb=0.0
      tot_nhb=0
      
      hb(:)=0.0
      nhb_exist(:)=0 
      nhb_start(:)=0 
      !=============
      !The main loop
      !=============     
      kLOOP: do k=1,nwat
        qj=0
        nqj=0 ! The number of bonded times for k-th form of quasi-HB 
        m1=ndx_1(k)
        m2=ndx_2(k)
        ndx_3_list=h_ndx_list(ndx_1(k),ndx_2(k),pos_filename,nat,boxsize)
        WRITE(*,*) "The ",k,"-th pair: ndx_of H (1st,2nd,3rd,4th):",& 
            ndx_3_list(1),ndx_3_list(2),ndx_3_list(3),ndx_3_list(4)
        ! Calculate h(j)
        ! A LOOP on ndx_3_list
        TIME: do jj=1,nmo
          h(jj)=0.0
          hb_exist(jj)=.False.
          hb_start(jj)=.False.

          ! Check if the pairs are located in one of the interfaces 
          index_mol1=grid_index(atom_info(m1,jj)%coord(1),&
              atom_info(m1,jj)%coord(2),divx,divy,nb_divx) 
          WRITE(*,*) "index_mol1=",index_mol1
          index_mol2=grid_index(atom_info(m2,jj)%coord(1),&
            atom_info(m2,jj)%coord(2),divx,divy,nb_divx) 

          !For surf 1
          condition1=pair_in_surf1(surf_info(index_mol1,jj)%coord(1),&
              atom_info(m1,jj)%coord(3),&
              surf_info(index_mol2,jj)%coord(1),&
              atom_info(m2,jj)%coord(3),thickness ) 

          !For surf 2 
          condition2=pair_in_surf2(surf_info(index_mol1,jj)%coord(2),&
              atom_info(m1,jj)%coord(3),&
              surf_info(index_mol2,jj)%coord(2),&
              atom_info(m2,jj)%coord(3),thickness ) 
          WRITE(*,*) condition1,condition2 

          !This condition is the additional condition for the establishment 
          ! of interface hydrogen bonds,which is the core of this method. 
          IF (condition1 .OR. condition2) THEN

            HYDROGEN: DO j=1,4
              m3=ndx_3_list(j)
              r1=(/atom_info(m1,jj)%coord(1),atom_info(m1,jj)%coord(2),&
                     atom_info(m1,jj)%coord(3) /)
              r2=(/atom_info(m2,jj)%coord(1),atom_info(m2,jj)%coord(2),&
                     atom_info(m2,jj)%coord(3) /)
              r3=(/atom_info(m3,jj)%coord(1),atom_info(m3,jj)%coord(2),&
                     atom_info(m3,jj)%coord(3) /)
              r21=dist2(r1,r2,boxsize) 
              if (criterion==1) then
                r23=dist2(r3,r2,boxsize) 
                r13=dist2(r3,r1,boxsize) 
                pm=pmADH(r1,r2,r3,boxsize)  ! if H is bound to O2
                pm_=pmADH(r2,r1,r3,boxsize) ! if H is bound to O1
                norm_rr=sqrt(r21*r23)
                cosphi=pm/norm_rr
                cosphi_=pm_/norm_rr
                if ((r21 .lt. rooc ).and. ( (cosphi .gt. cosPhiC123) .or. &
                    (cosphi_ .gt. cosPhiC123) )                      &
                   ) then
                  WRITE(*,*) "# of HB along time: ",qj+1    
                  h(jj)=1.0
                  hb_exist(jj)=.True.
                  !SPECIAL FOR S_HB(t) 
                  IF (jj .eq. 1) THEN
                    hb_start(jj)=.True.
                    qj=qj+h(jj)
                    nqj=nqj+1
                    EXIT ! if we know that two pair of molecule is bonded at frame jj,then we go to check the next frame (jj+1)
                  ELSE
                    IF(h(jj-1)<h_min)THEN
                      hb_start(jj)=.True. !If a bond is broken,recount a hydrogen bond
                      qj=qj+h(jj)
                      nqj=nqj+1 
                      EXIT ! if we know that two pair of molecule is bonded at frame jj,then we go to check the next frame (jj+1)
                    ENDIF
                  ENDIF
                endif
              elseif (criterion==2) then
                !Follow the second cirterion of HB.
                r31=dist2(r1,r3,boxsize) 
                r32=dist2(r2,r3,boxsize) 
                pm=pmAHD(r1,r2,r3,boxsize)
                cosphi=pm/(sqrt(r31*r32))

                !Follow the scond criterion of HB.
                !-0.5 comes from cos(2pi/3)
                if (r21 .lt. rooc .and. cosphi .lt. cosPhiC132) then 
                  h(jj)=1.0 
                  hb_exist(jj)=.True.
                  !SPECIAL FOR S_HB(t) 
                  IF (jj .eq. 1) THEN
                    hb_start(jj)=.True.
                    qj=qj+h(jj)
                    nqj=nqj+1
                    EXIT ! if we know that two pair of molecule is bonded at frame jj,then we go to check the next frame (jj+1)
                  ELSE
                    IF(h(jj-1)<h_min)THEN
                      hb_start(jj)=.True. !If a bond is broken,recount a hydrogen bond
                      qj=qj+h(jj)
                      nqj=nqj+1 
                      EXIT ! if we know that two pair of molecule is bonded at frame jj,then we go to check the next frame (jj+1)
                    ENDIF
                  ENDIF
                endif
              endif
            END DO HYDROGEN
          ENDIF
        enddo TIME 
        hb(k)=qj 
        nhb_exist(k)=nqj
        nhb_start(k)=nqj
        tot_hb=tot_hb + qj
        tot_nhb=tot_nhb+nhb_start(k) ! start
        !==================================
        !Calcualte the correlation function
        !==================================
        if (hb(k)>hb_min) then
          do mt=0,nmo-1    ! time interval
            scalar_hh=0.d0
            !do j=1,nmo-mt-1
            do j=1,nmo-mt
            hh=0.0
              !==================================================================== 
              !The definition/calculation of H(t): hh(j+mt)=h(j)*h(j+1)*...*h(j+mt)
              IF (hb_start(j))THEN
                DO m=0,mt
                  IF (h(j+m)>hb_min) THEN
                    hh=1.0
                  ELSE
                    hh=0.0
                    EXIT
                  ENDIF
                ENDDO
              ENDIF  
              !End of definition/calculation of H(t).
              !==================================================================== 
            scalar_hh=scalar_hh+h(j)*hh  
            enddo
          corr_hh(mt+1)=corr_hh(mt+1)+scalar_hh !sum_S_k(t)
          enddo
        endif
      enddo kLOOP   
      deallocate(hb_exist,nhb_exist)
      hb_per_frame=tot_hb/nmo
      WRITE(*,*) "Total H-bonds exists in history: ",tot_hb
      ave_h=hb_per_frame/nwat 
      WRITE(*,*) "hb per frame:",hb_per_frame
      WRITE(*,*) "<h>=",ave_h
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
      !Normalization of S_HB(t) step1
      !==============================
      do mt=0,nmo-1! time interval
          corr_hh(mt+1)=corr_hh(mt+1)/(nmo-mt)
      enddo
      corr_hh=corr_hh/nwat
      !Normalization step2
      corr_hh=corr_hh/ave_h
      deallocate(x,y,z,ndx_1,ndx_2)          
     !===================================
     !Write the correlation
     !S_HB(t) for the iterfacial HB (ihb)    
     !===================================
      char_thickness=nth(str(nint(thickness)),1)
      WRITE(char_criterion,'(i1)') criterion
      OPEN(10,file=trim(filename)//&
        '_'//donor//accepter//'_hbacf_hh_ihb_'//&
        char_thickness//"_"//char_criterion//'.dat')
        do i=1,nmo
            WRITE(10,*)(i-1)*delta_t,corr_hh(i)
        enddo
        WRITE(6,*)'written in '//trim(filename)//&
        '_'//donor//accepter//'_hbacf_hh_ihb_'//char_thickness//&
        "_"//char_criterion//'.dat'
      close(10)
     !=====================
     !Write the averate lifetime
     !<tau_a>
     !=====================
      tau_a=0.0
      OPEN(10,file=trim(filename)//&
        '_'//donor//accepter//'_hbacf_tau_a_shb_'//&
        char_thickness//&
        "_"//char_criterion//'.dat')
        do i=1,nmo
            tau_a=tau_a + delta_t * corr_hh(i)
        enddo
        WRITE(10,*) "HB lifetime (tau_a) (unit: ps): ", tau_a 
        WRITE(6,*)'written in '//trim(filename)//&
          '_'//donor//accepter//'_hbacf_tau_a_shb_'//char_thickness//&
          "_"//char_criterion//'.dat'
      close(10)
     !=====================
     !Write the correlation
     !ln(S_HB(t))     
     !=====================
      OPEN(10,file=trim(filename)//&
        '_'//donor//accepter//'_hbacf_log_hh_ihb_'//&
        char_thickness//"_"//char_criterion//'.dat')
        do i=1,nmo
            WRITE(10,*)(i-1)*delta_t,log(corr_hh(i))
        enddo
        WRITE(6,*)'written in '//trim(filename)//&
        '_'//donor//accepter//'_hbacf_log_hh_ihb_'//char_thickness//&
        "_"//char_criterion//'.dat'
      close(10)
     !===========
     ! Print <h>      
     !===========     
      OPEN(10,file=trim(filename)//&
        '_'//donor//accepter//'_ave_h_by_S_ihb_'//&
        char_thickness//"_"//char_criterion//'.dat')
        WRITE(10,*) 'Ave. No. bonds:',hb_per_frame
        WRITE(10,*) '<h>:',ave_h
        WRITE(6,*)'written in '//trim(filename)//&
        '_'//donor//accepter//'_ave_h_by_S_ihb_'//char_thickness//&
        "_"//char_criterion//'.dat'
      close(10)
      deallocate (h,corr_hh,hb)
      END SUBROUTINE 
