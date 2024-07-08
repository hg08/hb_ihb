      SUBROUTINE ghbacf_interface_n_pbc_format2(boxsize,delta_t0, &
          filename,pos_filename,list_filename,n_samples,nat,ns, &
          criterion,atom_info,n_grid,divx,divy,divz,nb_divx,nb_divy,&
          nb_divz,thickness,surf_info)
      !20.d0-5-29: simplifying the definition of hb (without r13)
      !20.d0-9-18: In this function, we define interfacial HB. 
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
      use module_ihb, ONLY: get_total_number_of_lines, &
                       hydrogen_ndx_list, &
                       distance2, &
                       pm_adh, pm_ahd, & 
                       grid_index, pair_in_surf1, pair_in_surf2,&
                       str, nth, atom
      implicit none
      !==========
      !parameters
      !==========
      INTEGER,parameter :: rk=8 ! local 
      INTEGER, parameter :: d_len=3 ! for storing the length of the character which represents the thickness of the interface

      character(LEN=200),INTENT(INOUT) :: filename,pos_filename
      character(LEN=200),INTENT(IN) :: list_filename
      INTEGER,INTENT(IN) :: criterion
      INTEGER, INTENT(IN) :: nat ! number of atoms
      INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
      REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: boxsize
      TYPE(atom), DIMENSION(nat,n_samples), INTENT(IN) :: atom_info
      REAL(KIND=rk), INTENT(IN) :: thickness ! the thickness of the instantaneous interfaces
      INTEGER,INTENT(IN) :: nb_divx, nb_divy, nb_divz, n_grid 
      REAL(KIND=rk),INTENT(IN) :: divx, divy, divz
      REAL(KIND=rk),DIMENSION(2,n_grid,n_samples),INTENT(IN) :: &
          surf_info
    
      !Local variables
      REAL(KIND=rk),parameter :: rooc=12.25d0                 ! cutoff distance of rOO (3.5**2 )
      REAL(KIND=rk),parameter :: cosPhiC123=0.866d0           ! 1.732/2; phiC=pi/6.
      REAL(KIND=rk),parameter :: cosPhiC132=-0.5d0            ! -1./2; phiC132=2pi/3.
      !REAL(KIND=rk), parameter :: h_min=0.d5  ! condition for the existence of a h-bond for a step
      REAL(KIND=rk),parameter :: hb_min=0.5d0 ! condition for the existence of h-bond for a pair of water molecules
      REAL(KIND=rk) :: r13, cosphi, pm, cosphi_, pm_
      REAL(KIND=rk) :: r21, r31, r32, r23 ! For the second criterion of HB
      REAL(KIND=rk) ::qj, tot_hb, delta_t, delta_t0, hb_per_frame, ave_h
      REAL(KIND=rk), DIMENSION(3) :: r1, r2, r3 ! pbc 
      INTEGER :: m1,m2,m3,mt,nqj,tot_nhb,n_bonded_pairs,ns
      REAL(KIND=rk),ALLOCATABLE,DIMENSION (:)  :: h,h_d,hb,corr_n
      REAL(KIND=rk),ALLOCATABLE,DIMENSION (:,:) :: x,y,z
      INTEGER,ALLOCATABLE,DIMENSION(:) :: ndx_1,ndx_2,nhb_exist
      INTEGER,DIMENSION(4) :: ndx_3_list
      REAL(KIND=rk)  :: scalar 
      logical,ALLOCATABLE,DIMENSION (:)  :: hb_exist
      INTEGER  :: nmo  ! nmo is not necessary, we set nmo = n_samples, because we do not want to change too much
      INTEGER :: nwat ! number of water molecules
      INTEGER :: i,j,k,jj 
      CHARACTER(len=d_len) :: char_thickness ! for saving the thickness in the files' names
      INTEGER :: index_mol1, index_mol2
      LOGICAL :: condition1,condition2
      !==============
      !Initialization
      !==============
      ave_h =0.d0; scalar =0.d0
      pm =0.d0; cosphi =0.d0
      r21 = 0.d0; r23 = 0.d0
      r31 = 0.d0; r13= 0.d0; r32 = 0.d0
      hb_per_frame = 0.d0; tot_hb = 0.d0
      r1 = 0.d0; r2 = 0.d0; r3 = 0.d0
      nmo = n_samples; nwat=0 
      ndx_3_list=0
      index_mol1=0; index_mol2=0
      condition1=.FALSE.
      condition2=.FALSE.

      !To obtain the total number of water pairs
      nwat=get_total_number_of_lines(list_filename)
      WRITE(*,*) 'ghbacf_c: # of water pairs (nwat) =', nwat
      ALLOCATE(ndx_1(nwat))          
      ALLOCATE(ndx_2(nwat))          
      !============================
      !read data from the list file
      !============================
      OPEN(10,file=list_filename, STATUS='OLD', ACTION='READ')     
      DO k=1,nwat
          read(10,*)ndx_1(k),ndx_2(k)
      ENDDO 
      CLOSE(10)
      !============================

      delta_t=REAL(ns,rk)*delta_t0  ! unit: ps
      WRITE(*,*) "New total steps (nmo):", nmo
      ALLOCATE(x(nat,nmo))
      ALLOCATE(y(nat,nmo))
      ALLOCATE(z(nat,nmo))
      ALLOCATE(h(nmo))
      ALLOCATE(h_d(nmo))
      ALLOCATE(hb(nwat)) ! Average H-bonded population 
      ALLOCATE(nhb_exist(nwat))
      !============================
      !read in surf trajectory file 
      !============================
      OPEN(20,file='surf_traj.dat')
           
      !==============================      
      ! Calculate correlation n_HB(t)  
      !==============================      
      ALLOCATE(corr_n(nmo))
      ALLOCATE(hb_exist(nmo))
      corr_n(:)=0.d0      
      tot_hb=0.d0
      tot_nhb=0
      h=0.d0; h_d=0.d0
      hb(:)=0.d0
      nhb_exist(:)=0 
      !=============
      !The main loop
      !=============      
      kLOOP: DO k=1,nwat
        qj=0
        nqj=0 ! The number of bonded times for k-th form of quasi-HB 
        m1=ndx_1(k)
        m2=ndx_2(k)
        ndx_3_list= & 
            hydrogen_ndx_list(ndx_1(k),ndx_2(k),pos_filename,nat,boxsize)
        !WRITE(*,*) "The ",k,"-th pair: ndx_of H (1st,2nd,3rd,4th):",& 
        !    ndx_3_list(1), ndx_3_list(2), ndx_3_list(3), ndx_3_list(4)
        ! Calculate h_d(j)
        TIME: DO jj =1, nmo
          h(jj)=0.d0
          hb_exist(jj)=.False.

          ! Check if the pairs are located in one of the interfaces 
          index_mol1 =grid_index(atom_info(m1,jj)%coord(1), &
              atom_info(m1,jj)%coord(2),divx,divy,nb_divx,nb_divy) 
          !WRITE(*,*) "index_mol1 = ",index_mol1
          index_mol2 =grid_index(atom_info(m2,jj)%coord(1), &
              atom_info(m2,jj)%coord(2),divx,divy,nb_divx,nb_divy) 

          !For surf 1
          condition1 = pair_in_surf1(surf_info(1,index_mol1,jj),&
              atom_info(m1,jj)%coord(3), &
              surf_info(1,index_mol2,jj), &
              atom_info(m2,jj)%coord(3),thickness ) 

          !For surf 2 
          condition2 = pair_in_surf2(surf_info(2,index_mol1,jj),&
              atom_info(m1,jj)%coord(3), &
              surf_info(2,index_mol2,jj), &
              atom_info(m2,jj)%coord(3),thickness ) 
          !WRITE(*,*) condition1, condition2 

          !This condition is the additional condition for the establishment 
          ! of interface hydrogen bonds, which is the core of this method. 
          IF (condition1 .OR. condition2) THEN

              ! A LOOP on ndx_3_list
              HYDROGEN: DO j=1,4
                  !Try each Hydrogen bond, there are totally 4 hydrogen atoms, for a water-water pair
                  m3=ndx_3_list(j)
                  r1 = (/atom_info(m1,jj)%coord(1), &
                         atom_info(m1,jj)%coord(2), &
                         atom_info(m1,jj)%coord(3)/)
                  r2 = (/atom_info(m2,jj)%coord(1), &
                         atom_info(m2,jj)%coord(2), &
                         atom_info(m2,jj)%coord(3)/)
                  r3 = (/atom_info(m3,jj)%coord(1), &
                         atom_info(m3,jj)%coord(2), &
                         atom_info(m3,jj)%coord(3)/)
                  r21 = distance2(r1, r2, boxsize)

                  IF (criterion==1) THEN
                      !Follow the first cirterion of HB.
                      r23 = distance2(r3, r2, boxsize) 
                      r13 = distance2(r3, r1, boxsize) 
                      pm = pm_adh(r1, r2, r3, boxsize)  ! IF H is bound to O2
                      pm_ = pm_adh(r2, r1, r3, boxsize) ! If H is bound to O1
                      cosphi = pm/(sqrt(r21*r23))
                      cosphi_= pm_/(sqrt(r21*r23))
                      IF (r21 .lt. rooc ) THEN
                          h_d(jj) = 1.0d0 
                      ENDIF
                      IF ( (r21 .lt. rooc).and.((cosphi.gt.cosPhiC123) .or. &
                          (cosphi_ .gt. cosPhiC123) )        &
                         ) THEN
                          h(jj)=1.0d0 
                          hb_exist(jj) = .True.
                          qj=qj+h(jj) ! To calculate ave population of HB over all starting points for one pair of water                           
                          nqj=nqj+1
                          EXIT ! If we know that two pair of molecule is bonded at frame jj, then we go to check the next frame (jj+1)
                      ENDIF 
                  ELSEIF (criterion==2) THEN
                      !Follow the second cirterion of HB.
                      r31 = distance2(r1,r3,boxsize) 
                      r32 = distance2(r2,r3, boxsize) 
                      pm= pm_ahd(r1,r2,r3,boxsize)
                      cosphi= pm/(sqrt(r31*r32))
                      IF (r21 .lt. rooc ) THEN
                          h_d(jj)=1.0d0 
                      ENDIF
                      IF ((r21 .lt. rooc ).and.(cosphi .lt. cosPhiC132) &
                         ) THEN
                          h(jj)=1.0d0 
                          hb_exist(jj) = .True.
                          qj=qj+h(jj) ! To calculate ave population of HB over all starting points for one pair of water                           
                          nqj=nqj+1
                          EXIT ! If we know that two pair of molecule is bonded at frame jj, then we exit and go to check the next frame (jj+1)
                      ENDIF 
                   ENDIF
               END DO HYDROGEN
           ENDIF
        ENDDO   TIME 
        hb(k)=qj 
        nhb_exist(k)=nqj
        tot_hb=tot_hb + qj
        tot_nhb=tot_nhb+nhb_exist(k)
        !==========================================
        !Calcualte the correlation function n_HB(t)
        !==========================================
        IF (hb(k)>hb_min) THEN
            DO mt=0,nmo-1    ! time interval
                scalar=0.0d0
                !DO j=1,nmo-mt-1
                DO j=1,nmo-mt
                    scalar=scalar+h(j)*(1-h(j+mt))*h_d(j+mt)  
                ENDDO  
                corr_n(mt+1)=corr_n(mt+1)+scalar !sum_C_k(t)
            ENDDO 
        ENDIF
      ENDDO  kLOOP   
      DEALLOCATE(hb_exist,nhb_exist,h_d)
      hb_per_frame = tot_hb/REAL(nmo,rk)
      !WRITE(*,*) "Total H-bonds exists in history: ", tot_hb
      ave_h = hb_per_frame/REAL(nwat,rk) 
      !WRITE(*,*) "hb per frame:", hb_per_frame
      !WRITE(*,*) "<h> =", ave_h
      !=========================================
      !Calculate the number of ever bonded pairs
      !=========================================
      n_bonded_pairs=0 
      DO k=1,nwat
          IF (hb(k)>hb_min) THEN
              n_bonded_pairs=n_bonded_pairs+1      
          ENDIF
      ENDDO 
      !==============================
      !Normalization of n_HB(t) step1
      !==============================
      DO mt=0,nmo-1! time interval
          corr_n(mt+1)=corr_n(mt+1)/REAL(nmo-mt,rk)
      ENDDO 
      corr_n = corr_n /REAL(nwat,rk)
      !Normalization step2
      corr_n = corr_n/ave_h
      DEALLOCATE(x,y,z,ndx_1,ndx_2)          
     !====================================
     !WRITE the correlation
     !n_HB(t) for the interfacial HB (ihb)     
     !====================================
      !char_thickness = nth(str(nint(thickness)),d_len)
      char_thickness = nth(str(thickness),d_len)
      OPEN(10,file=trim(filename)//'_wat_pair_hbacf_n_ihb_'//&
        char_thickness//'.dat')
        DO i=1,nmo
            WRITE(10,*) REAL(i-1,rk) * delta_t,corr_n(i)
        ENDDO 
        WRITE(6,*)'written in '//trim(filename)//&
                  '_wat_pair_hbacf_n_ihb_'//char_thickness//'.dat'
      CLOSE(10)
     !=====================
     !WRITE the correlation
     !ln(n_HB(t))     
     !=====================
      OPEN(10,file=trim(filename)//'_wat_pair_hbacf_log_n_ihb_'//&
        char_thickness//'.dat')
        DO i=1,nmo
            WRITE(10,*) REAL(i-1,rk)*delta_t,log(corr_n(i))
        ENDDO 
        WRITE(6,*)'written in '//trim(filename)//&
                  '_wat_pair_hbacf_log_n_ibh_'//char_thickness//'.dat'
      CLOSE(10)
     !===========
     ! Print <h>      
     !===========      
      OPEN(10,file=trim(filename)//'_wat_pair_ave_h_by_n_ihb_'//&
        char_thickness//'.dat')
        WRITE(10,*) 'Ave. No. bonds:', hb_per_frame
        WRITE(10,*) '<h>:',ave_h
        WRITE(6,*)'written in '//trim(filename)//&
                  '_wat_pair_ave_h_by_n_ihb_'//char_thickness//'.dat'
      CLOSE(10)
      DEALLOCATE (h,corr_n,hb)
      END SUBROUTINE 
