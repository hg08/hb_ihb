      SUBROUTINE ghbacf_interface_n_pbc_format2(boxsize, delta_t0, &
          filename, pos_filename, list_filename, n_samples, nat, ns, &
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
      !5)A local PARAMETER thickness is needed to define the thickness of 
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
                       grid_index, pair_in_surf1, pair_in_surf2, &
                       str, nth
      USE atom_module
      USE parameter_shared ! Including ndx_O, ndx_H, ...
      USE traj_format2
      USE surf_traj
      USE surf_module, ONLY: surf_info

      IMPLICIT NONE
      !==========
      !PARAMETERS
      !==========
      INTEGER,PARAMETER :: rk=8 ! local 
      INTEGER, PARAMETER :: d_len=1 ! for storing the length of the character which represents the thickness of the interface
      CHARACTER(LEN=200),INTENT(INOUT) :: filename, pos_filename
      CHARACTER(LEN=200),INTENT(IN) :: list_filename
      INTEGER,INTENT(IN) :: criterion
      INTEGER,INTENT(IN) :: nat ! number of atoms
      INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
      REAL(KIND=rk),DIMENSION(3), INTENT(IN) :: boxsize
      REAL(KIND=rk),INTENT(IN) :: delta_t0
    
      !Local variables
      REAL,PARAMETER :: rooc=12.25 ! cutoff distance of rOO (3.5**2 )
      REAL,PARAMETER :: cosPhiC123=0.866 ! 1.732/2; phiC=pi/6.
      REAL,PARAMETER :: cosPhiC132=-0.5 ! -1./2; phiC132=2pi/3.
      REAL(KIND=rk),PARAMETER :: max_time_for_corr = 12.0 ! Unit: ps. 

      REAL(KIND=rk),PARAMETER :: h_min=0.5 ! condition for the existence of a h-bond for a step
      REAL(KIND=rk),PARAMETER :: hb_min=0.5 ! condition for the existence of h-bond for a pair of water molecules
      REAL(KIND=rk) :: r13, cosphi, pm, cosphi_, pm_, norm_rr
      REAL(KIND=rk) :: r21, r31, r32, r23 ! For the second criterion of HB
      REAL(KIND=rk) :: qj, tot_hb, delta_t, hb_per_frame, ave_h
      REAL(KIND=rk),DIMENSION(3) :: r1, r2, r3 ! pbc 
      INTEGER :: m1, m2, m3, mt, nqj, tot_nhb, n_bonded_pairs, ns
      INTEGER :: idx_O1, idx_O2 ! Indices of O1, O2 in all O atoms
      REAL(KIND=rk),ALLOCATABLE,DIMENSION (:) :: h, hb, corr_n
      REAL(KIND=rk),ALLOCATABLE,DIMENSION (:) :: h_d
      REAL,ALLOCATABLE,DIMENSION (:,:) :: x, y, z
      INTEGER,ALLOCATABLE,DIMENSION(:) :: nhb_exist
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ndx_1, ndx_2
      INTEGER :: nmo_effective, start_step, num_start_points
      INTEGER,DIMENSION(4) :: ndx_3_list
      REAL(KIND=rk) :: scalar, tmp 
      LOGICAL,ALLOCATABLE,DIMENSION (:) :: hb_exist
      INTEGER :: nmo ! nmo is not necessary, we set nmo = n_samples, because we DO not want to change too much
      INTEGER :: nwat ! number of water molecules
      INTEGER :: i, j, k, jj 
      CHARACTER(LEN=d_len) :: char_thickness ! for saving the thickness in the files' names
      INTEGER :: index_mol1, index_mol2
      LOGICAL :: condition1, condition2
      !==============
      !Initialization
      !==============
      ave_h = 0.0; scalar = 0.0
      pm = 0.0; cosphi = 0.0
      r21 = 0.0; r23 = 0.0
      r31 = 0.0; r13 = 0.0; r32 = 0.0
      hb_per_frame = 0.0; tot_hb = 0.0
      r1 = 0.0; r2 = 0.0; r3 = 0.0
      nmo = n_samples; nwat = 0 
      ndx_3_list = 0
      index_mol1 = 0; index_mol2 = 0
      condition1 = .FALSE.
      condition2 = .FALSE.
      norm_rr = 0.0 ! a temporary variable
      tmp = 0.0 ! a temporay variable 
      char_thickness = ''
      start_step = 1
      nmo_effective = 0

      !To obtain the total number of water pairs
      nwat = get_nwat(list_filename)
      ALLOCATE(ndx_1(nwat,2))          
      ALLOCATE(ndx_2(nwat,2))          
      !============================
      !read data from the list file
      !============================
      OPEN(10,file=list_filename)     
      DO k = 1, nwat
          read(10,*)ndx_1(k,1), ndx_1(k,2), ndx_2(k,1), ndx_2(k,2)
      ENDDO
      CLOSE(10)
      !============================

      delta_t = ns * delta_t0  ! unit: ps
      nmo_effective = nint(max_time_for_corr/delta_t) + 1
      start_step = nint((nmo_effective-1)/10.0) ! Start step of sliding window. Over-using rate is 1 - 1/5 = 4/5
      num_start_points = (nmo-nmo_effective-1)/start_step + 1
      WRITE(*,*) "New time step (delta_t):", delta_t
      WRITE(*,*) "New total steps (nmo):", nmo
      ALLOCATE(x(nat,nmo))
      ALLOCATE(y(nat,nmo))
      ALLOCATE(z(nat,nmo))
      ALLOCATE(h(nmo))
      ALLOCATE(h_d(nmo))
      ALLOCATE(hb(nwat))    ! Average H-bonded population 
      ALLOCATE(nhb_exist(nwat))
      !============================
      !read in surf trajectory file 
      !============================
      !OPEN(20,file='surf_traj.dat')
           
      !==============================      
      ! Calculate correlation n_HB(t)  
      !==============================      
      ALLOCATE(corr_n(nmo))
      ALLOCATE(hb_exist(nmo))
      corr_n(:) = 0.0      
      tot_hb = 0.0
      tot_nhb = 0
      h = 0.0
      h_d = 0.0
      hb(:) = 0.0
      nhb_exist(:) = 0 
      !=============
      !The main loop
      !=============      
      kLOOP: DO k = 1, nwat ! Loop all O-O pairs
        qj = 0
        nqj = 0 ! The number of bonded times for k-th form of quasi-HB 
        idx_O1 = ndx_1(k,1) ! Index of O1 in all O atoms
        idx_O2 = ndx_2(k,1) ! Index of O2 in all O atoms
        m1 = ndx_1(k,2) ! Index of O1 in all atoms
        m2 = ndx_2(k,2) ! Index of O2 in all atoms
        ndx_3_list = (/ndx_H(2*idx_O1-1,2), ndx_H(2*idx_O1,2), ndx_H(2*idx_O2-1,2), ndx_H(2*idx_O2,2)/)
        ! Calculate h_d(j)
        TIME: DO jj = 1, nmo
          h(jj) = 0.0
          hb_exist(jj) = .False.

          ! Check IF the pairs are located in one of the interfaces 
          index_mol1 = grid_index(atom_info(m1,jj)%coord(1), &
              atom_info(m1,jj)%coord(2),divx,divy,nb_divx) 
          index_mol2 = grid_index(atom_info(m2,jj)%coord(1), &
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

          !The following condition establish interface hydrogen bonds,
          !which is the core of this method.
          IF (condition1 .OR. condition2) THEN
              r1 = (/atom_info(m1,jj)%coord(1),atom_info(m1,jj)%coord(2),&
                     atom_info(m1,jj)%coord(3) /)
              r2 = (/atom_info(m2,jj)%coord(1),atom_info(m2,jj)%coord(2),&
                     atom_info(m2,jj)%coord(3) /)
              HYDROGEN: DO j = 1, 4
                  !Try each Hydrogen bond, there are totally 4 hydrogen atoms, for a water-water pair
                  m3 = ndx_3_list(j)
                  r3 = (/atom_info(m3,jj)%coord(1),atom_info(m3,jj)%coord(2),&
                         atom_info(m3,jj)%coord(3) /)
                  r21 = dist2(r1, r2, boxsize)

                  IF (criterion==1) THEN
                      !Follow the first cirterion of HB.
                      r23 = dist2(r3, r2, boxsize) 
                      r13 = dist2(r3, r1, boxsize) 
                      pm = pmADH(r1,r2,r3,boxsize)  ! IF H is bound to O2
                      pm_ = pmADH(r2,r1,r3,boxsize) ! IF H is bound to O1
                      norm_rr = sqrt(r21*r23)
                      cosphi = pm/norm_rr
                      cosphi_= pm_/norm_rr
                      IF (r21 .lt. rooc ) THEN
                          h_d(jj) = 1.0 
                      ENDIF
                      IF ( (r21 .lt. rooc) .and. ( (cosphi .gt. cosPhiC123) .or. &
                          (cosphi_ .gt. cosPhiC123) )        &
                         ) THEN
                          h(jj) = 1.0 
                          hb_exist(jj) = .True.
                          qj = qj + h(jj) ! To calculate ave population of HB over all starting points for one pair of water                           
                          nqj = nqj + 1
                          EXIT ! IF we know that two pair of molecule is bonded at frame jj, THEN we go to check the next frame (jj+1)
                      ENDIF 
                  ELSEIF (criterion==2) THEN
                      !Follow the second cirterion of HB.
                      r31 = dist2(r1,r3,boxsize) 
                      r32 = dist2(r2,r3, boxsize) 
                      pm = pmAHD(r1,r2,r3,boxsize)
                      cosphi = pm/(sqrt(r31*r32))
                      IF (r21 .lt. rooc ) THEN
                          h_d(jj) = 1.0 
                      ENDIF
                      IF ((r21 .lt. rooc ).and.(cosphi .lt. cosPhiC132) &
                         ) THEN
                          h(jj) = 1.0 
                          hb_exist(jj) = .True.
                          qj = qj + h(jj) ! To calculate ave population of HB over all starting points for one pair of water                           
                          nqj = nqj + 1
                          EXIT ! IF we know that two pair of molecule is bonded at frame jj, THEN we exit and go to check the next frame (jj+1)
                      ENDIF 
                   ENDIF
               END DO HYDROGEN
           ENDIF
        ENDDO TIME 
        hb(k) = qj 
        nhb_exist(k) = nqj
        tot_hb = tot_hb + qj
        tot_nhb = tot_nhb + nhb_exist(k)
        !==========================================
        !Calcualte the correlation function n_HB(t)
        !==========================================
        IF (hb(k)>hb_min) THEN
            DO mt =0, nmo_effective-1   ! time interval
                scalar = 0.0
                DO j = 1, nmo-nmo_effective, start_step ! How many steps? (nmo-nmo_effective-1)/start_step + 1 
                    tmp = h(j) * (1-h(j+mt)) * h_d(j+mt)
                    scalar = scalar + tmp  
                ENDDO
                corr_n(mt+1) = corr_n(mt+1) + scalar !sum_n_k(t)
            ENDDO
        ENDIF
      ENDDO kLOOP   
      deALLOCATE(hb_exist,nhb_exist,h_d)
      hb_per_frame = tot_hb/REAL(nmo,rk)
      ave_h = hb_per_frame/REAL(nwat,rk) 
      !=========================================
      !Calculate the number of ever bonded pairs
      !=========================================
      n_bonded_pairs = 0 
      DO k = 1, nwat
          IF (hb(k)>hb_min) THEN
              n_bonded_pairs = n_bonded_pairs+1      
          ENDIF
      ENDDO
      !==============================
      !Normalization of n_HB(t) step1
      !==============================
      !DO mt = 0, nmo-1! time interval
      !    corr_n(mt+1) = corr_n(mt+1)/num_start_points
      !ENDDO
      corr_n = corr_n/(num_start_points * nwat)
      !Normalization step2
      corr_n = corr_n / ave_h
      deALLOCATE(x,y,z,ndx_1,ndx_2)          
     !====================================
     !Write the correlation
     !n_HB(t) for the interfacial HB (ihb)     
     !====================================
      char_thickness = nth(str(nint(thickness)),d_len)
      OPEN(10,file=trim(filename)//'_wat_pair_hbacf_n_ihb_'//&
        char_thickness//'.dat')
        DO i = 1, nmo_effective 
            WRITE(10,*) REAL(i-1)*delta_t, corr_n(i)
        ENDDO
        WRITE(6,*)'written in '//trim(filename)//&
                  '_wat_pair_hbacf_n_ihb_'//char_thickness//'.dat'
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
      deALLOCATE (h, corr_n, hb)
      END SUBROUTINE 
