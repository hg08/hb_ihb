      SUBROUTINE ghbacf_interface_nf_pbc_format2(boxsize, delta_t0, &
          filename, pos_filename, list_filename,list_filename_H, n_samples, nat, ns, &
          criterion)
      !========================================================================
      !Purpose: to obtain C_feeOH(t):  <nf(0)nf(t)>/<nf> for interfacial free OH
      !1)the correlation function is the average over OH goups
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
      ! Hydrogen bond definition (1 or 2)
      ! surf_filename
      !=========================================================
      !Modules
      use tools, ONLY: get_nline => get_total_numer_of_lines, &
                       h_ndx_list => hydrogen_ndx_list, &
                       dist2 => distance2, &
                       pmADH, pmAHD, & 
                       grid_index, pair_in_surf1, pair_in_surf2, &
                       atom_in_surf1, atom_in_surf2, &
                       str, nth
      USE atom_module
      USE water_molecule_types
      USE parameter_shared ! Including ndx_O, ndx_H, ...
      USE traj_format2
      USE surf_traj
      USE surf_module, ONLY: surf_info
 
      implicit none
      
      !==========
      !PARAMETERS
      !==========
      INTEGER,PARAMETER :: rk=8 ! local 
      INTEGER,PARAMETER :: d_len=1 ! for storing the length of the character which represents the thickness of the interface
      CHARACTER(LEN=200),INTENT(INOUT) :: filename, pos_filename
      CHARACTER(LEN=200),INTENT(IN) :: list_filename
      CHARACTER(LEN=200),INTENT(IN) :: list_filename_H
      INTEGER,INTENT(IN) :: criterion
      INTEGER,INTENT(IN) :: nat ! number of atoms
      INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
      REAL(KIND=rk),DIMENSION(3),INTENT(IN) :: boxsize
      REAL(KIND=rk),INTENT(IN) :: delta_t0

      !Local variables
      REAL,PARAMETER :: rooc=12.25 ! cutoff distance of rOO (3.5**2 )
      REAL,PARAMETER :: cosPhiC123=0.866 ! 1.732/2; phiC=pi/6.
      REAL,PARAMETER :: cosPhiC123_freeOH= 0.6428 !  phiC= 50 degree. (Ref. J. Chem. Theory Comput. 2018, 14, 357−364)
      REAL,PARAMETER :: cosPhiC132=-0.5 ! -1./2; phiC132=2pi/3.
      REAL,PARAMETER :: cosPhiC132_freeOH = -0.342 ! ; phiC132= 110. (Ref. Tang, J. Chem. Theory Comput. 2018, 14, 357−364)
      REAL(KIND=rk),PARAMETER :: max_time_for_corr = 12.0 ! Unit: ps. 
      REAL(KIND=rk),PARAMETER :: h_min=0.5 ! condition for the existence of a h-bond for a step
      REAL(KIND=rk),PARAMETER :: hb_min=0.5 ! condition for the existence of h-bond for a pair of water molecules
      REAL(KIND=rk),PARAMETER :: nfb_min=0.5 ! condition for the existence of free OH for a  OH in water molecule
      REAL(KIND=rk) :: r13, cosphi, pm, cosphi_, pm_, norm_rr
      REAL(KIND=rk) :: r21, r31, r32, r23 ! For the second criterion of HB
      REAL(KIND=rk) :: qj, tot_hb, delta_t, hb_per_frame, ave_h
      REAL(KIND=rk) :: freeoh, tot_nfb, nfb_per_frame, ave_nf 
      REAL(KIND=rk),DIMENSION(3) :: r1, r2, r3 ! pbc 
      INTEGER :: m1, m2, m3, mt, nqj, tot_nhb, n_bonded_pairs, ns
      INTEGER ::                 nfreeoh, tot_nfreeoh, n_freeoh
      INTEGER :: m1_neighbor, n_neighbor
      INTEGER :: idx_O1, idx_O2 ! Indices of O1, O2 in all O atoms
      INTEGER :: idx_H ! Indices of H in all O atoms
      REAL(KIND=rk),ALLOCATABLE,DIMENSION (:) :: h, hb, corr_h 
      REAL(KIND=rk),ALLOCATABLE,DIMENSION (:) :: nf, nfb, corr_nf
      REAL,ALLOCATABLE,DIMENSION (:,:) :: x, y, z
      INTEGER,ALLOCATABLE,DIMENSION(:) :: nhb_exist
      INTEGER,ALLOCATABLE,DIMENSION(:) :: nfreeoh_exist
      INTEGER,ALLOCATABLE,DIMENSION(:) :: nf_exist
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ndx_1, ndx_2
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ndx_H1, ndx_O1
      INTEGER,DIMENSION(4) :: ndx_3_list
      REAL(KIND=rk) :: scalar, tmp 
      LOGICAL,ALLOCATABLE,DIMENSION (:) :: hb_exist
      LOGICAL,ALLOCATABLE,DIMENSION (:) :: freeoh_exist
      LOGICAL,ALLOCATABLE,DIMENSION (:) :: is_free
      INTEGER :: nmo ! nmo is not necessary, we set nmo = n_samples, because we DO not want to change too much
      INTEGER :: nmo_effective, start_step, num_start_points 
      INTEGER :: nwat ! number of water molecules
      INTEGER :: n_H ! number of OH groups; or num of H atoms
      INTEGER :: n_O ! number of O atoms
      INTEGER :: i, j, k, jj 
      INTEGER :: k_H, k_O, k_O1, k_O2
      CHARACTER(LEN=d_len) :: char_thickness ! for saving the thickness in the files' names
      INTEGER :: index_mol1, index_mol2
      INTEGER :: index_mol
      LOGICAL :: condition1, condition2
      INTEGER :: tmp_index 

      !==============
      !Initialization
      !==============
      ave_h = 0.0; scalar = 0.0
      pm = 0.0; cosphi = 0.0
      r21 = 0.0; r23 = 0.0
      r31 = 0.0; r13 = 0.0; r32 = 0.0
      hb_per_frame = 0.0; tot_hb = 0.0
      nfb_per_frame = 0.0; tot_nfb = 0.0
      r1 = 0.0; r2 = 0.0; r3 = 0.0
      nmo = n_samples; nwat = 0 
      ndx_3_list = 0
      index_mol = 0
      index_mol1 = 0; index_mol2 = 0
      condition1 = .FALSE.
      condition2 = .FALSE.
      norm_rr = 0.0 ! a temporary variable
      tmp = 0.0 ! a temporay variable 
      char_thickness = ''
      start_step = 1
      nmo_effective = 0
      k_O = 0
      tmp_index = 0
       
      !To obtain the total number of water pairs
      nwat=get_nline(list_filename)
      n_H=get_nline(list_filename_H)
      n_O = n_H/2
      allocate(ndx_1(nwat,2))          
      allocate(ndx_2(nwat,2))          
      allocate(ndx_H1(n_H,2)) ! Indices of H in OH groups, ndx_H1(k_H, 1) is self-indices, ndx_H1(k_H, 2) is total-indices          
      allocate(ndx_O1(n_O, 2) ) ! Indices of O in OH groups, ndx_O1(k_O, 1) is self-indices, ndx_O1(k_O, 2) is total-indices          
      !============================
      !read data from the list file
      !============================
      OPEN(10, file=list_filename)     
      DO k = 1, nwat
          read(10,*)ndx_1(k,1), ndx_1(k,2), ndx_2(k,1), ndx_2(k,2)
      ENDDO
      close(10)
      OPEN(10, file=list_filename_H)     
      DO k_H = 1, n_H
          read(10,*)ndx_H1(k_H,1), ndx_H1(k_H,2)
      ENDDO
      close(10)
      !============================

      delta_t = ns * delta_t0 ! unit: ps
      nmo_effective = nint(max_time_for_corr/delta_t) + 1 
      start_step = nint((nmo_effective-1)/10.0) ! Start step of sliding window. Over-using rate is 1 - 1/5 = 4/5
      num_start_points = (nmo-nmo_effective-1)/start_step + 1 
      WRITE(*,*) "New total steps (nmo):", nmo
      allocate(x(nat,nmo))
      allocate(y(nat,nmo))
      allocate(z(nat,nmo))
      allocate(h(nmo))
      allocate(nf(nmo))
      allocate(hb(nwat)) ! Average H-bonded population 
      allocate(nhb_exist(nwat))
      allocate(nfreeoh_exist(n_H))
     !====================================
     ! Calculate <h(0)h(t)>/<h>  
     ! Notice here <> is not average over
     ! different pairs of water molecules,
     ! and over all starting time points i
     ! with h(i)=1.
     !====================================      
      allocate(corr_h(nmo))
      allocate(corr_nf(nmo))
      allocate(hb_exist(nmo))
      allocate(freeoh_exist(nmo))
      ! loop
      corr_h(:) = 0.0      
      tot_hb = 0.0
      tot_nhb = 0
      h(:) = 0.0 
      hb(:) = 0.0
      nhb_exist(:) = 0 
      !A loop over TOTAL Time Step
      TTLOOP: DO jj = 1,nmo
      DO k_O1 = 1, n_O
          ! TODO: 知道一个O,需要先计算其在任意一个时刻jj的O邻居序号之列表.
          m1 = ndx_O1(k_O, 2) ! Host O1, total index
          DO k_O2 = 1 , n_O
              if (k_O2 /= k_O1) then
                  m2 = ndx_O1(k_O2,2) ! For searching neighbors of Host O1
                  ! The total inices of H atoms related to O1 (Host O) and O2 (potential neighor O)
                  !ndx_3_list = (/ndx_H(2*k_O1-1,2), ndx_H(2*k_O1,2), ndx_H(2*k_O2-1,2), ndx_H(2*k_O2,2)/)
                  ! Check if O2 is a real neighbor of O1
                  r1 = (/atom_info(m1,jj)%coord(1),atom_info(m1,jj)%coord(2),&
                     atom_info(m1,jj)%coord(3) /)
                  r2 = (/atom_info(m2,jj)%coord(1),atom_info(m2,jj)%coord(2),&
                         atom_info(m2,jj)%coord(3) /)
                  r21 = dist2(r1, r2, boxsize)
                  if ((r21 .lt. rooc)) THEN
                     ! TODO: In the neighbor list of H atom, add the total indices of O2.
                     O_info(k_O1,jj)%num_oxygen_neighbors = O_info(k_O1,jj)%num_oxygen_neighbors + 1 ! increase the number of O neighbors
                     tmp_index = O_info(k_O1,jj)%num_oxygen_neighbors
                     O_info(O_info(k_O1,jj)%indices_oxygen_neighbors(tmp_index),jj)%atom_id = 3*(k_O2-1) + 1 ! define total index of O
                  endif
              endif          
          ENDDO
      ENDDO
      ENDDO TTLOOP

      tLOOP: DO jj = 1,nmo
          !===================================
          !The loop to find out free OH groups
          !===================================
          DO k_H = 1, n_H ! self index of H
              !ndx_O1(ceiling(ndx_H1(k_H,1)/2),1 ) ! Indices of self-order of Host O (O1)
              !ndx_O1(ceiling(ndx_H1(k_H,1)/2),2 ) ! Indices of total-order of Host O (O1)
              idx_H = H_info(k_H,jj)%atom_id ! total index of H
              m3 = idx_H 
              k_O1 =ceiling(k_H * 1.0 / 2) ! indices of self-indices of Host O
              m1 = H_info(k_H,jj)%host_id ! Indices of total index of the Host 
              qj = 0
              if (O_info(k_O1,jj)%num_oxygen_neighbors .le. 0) then
                  H_info(k_H,jj)%IsFree = .True.
              else
                  allocate(is_free(O_info(k_O1,jj)%num_oxygen_neighbors))
                  DO k_O2 = 1, O_info(k_O1,jj)%num_oxygen_neighbors 
                      m2 = O_info(k_O1,jj)%indices_oxygen_neighbors(k_O2) ! Indices of total index of the neighbor O2 
                      ! DONE: 一旦知此列表序列,对于任何一个氢原子,我们可以找到其宿主氧原子idx_O_host (H_info(:,:)%host_id)
                      ! 对于任意一对组合(idx_O_host, idx_O_neighbor[1]), (idx_O_host, idx_O_neighbor[2]), (idx_O_host, idx_O_neighbor[3]),...
                      !分别考察当前OH是否满足自由OH条件,随后根据结果综合判断该OH是否为自由OH.
                      ! so that we can Calculate nf(j) later.

                      r1 = (/atom_info(m1,jj)%coord(1),atom_info(m1,jj)%coord(2),&
                             atom_info(m1,jj)%coord(3) /)
                      r2 = (/atom_info(m2,jj)%coord(1),atom_info(m2,jj)%coord(2),&
                             atom_info(m2,jj)%coord(3) /)
                      r3 = (/atom_info(m3,jj)%coord(1),atom_info(m3,jj)%coord(2),&
                             atom_info(m3,jj)%coord(3) /)
                      r21 = dist2(r1, r2, boxsize) 
                      if (criterion == 1) THEN
                          r23 = dist2(r3, r2, boxsize) 
                          r13 = dist2(r3, r1, boxsize) 
                          pm_ = pmADH(r2,r1,r3,boxsize) ! if H is bound to O1
                          norm_rr = sqrt(r21*r23)
                          cosphi_= pm_/norm_rr
                          if (cosphi_  .lt. cosPhiC123_freeOH)  THEN
                              is_free(k_O2) = .True. ! The k_O2-th step for checking whether the k_H -th H is free.
                          endif
                      elseif (criterion == 2) THEN
                          !Follow the second cirterion of HB.
                          r31 = dist2(r1,r3,boxsize) 
                          r32 = dist2(r2,r3, boxsize) 
                          pm = pmAHD(r1,r2,r3,boxsize)
                          cosphi = pm/(sqrt(r31*r32))

                          !Follow the scond criterion of HB.
                          !-0.342 comes from cos(110 degree)
                          if (cosphi .gt. cosPhiC132_freeOH) THEN 
                              is_free(k_O2) = .True.
                          endif
                      endif
                      if ( all(is_free .eqv. .True.) .eqv. .True. ) then
                          H_info(k_H,jj)%IsFree = .True.
                      endif
                  ENDDO
              endif          
          ENDDO
      ENDDO tLOOP

      hLOOP: DO k_H = 1, n_H
        freeoh = 0
        nfreeoh = 0 ! The number of freeoh-ed times for k-th OH group
        k_O1 =ceiling(k_H * 1.0 / 2) ! indices of self-indices of Host O
        m1 = H_info(k_H,jj)%host_id ! Indices of total index of the Host (O1)
        DO jj = 1, nmo
          nf(jj) = 0.0
          freeoh_exist(jj) = .False.
          index_mol1 = grid_index(atom_info(m1,jj)%coord(1), &
              atom_info(m1,jj)%coord(2),divx,divy,nb_divx)  

          ! Check if O1 is in one of the two interfaces
          !For surf 1
          condition1 = atom_in_surf1(surf_info(index_mol,jj)%coord(1),&
                atom_info(m1,jj)%coord(3), thickness )
          !For surf 2 
          condition2 = atom_in_surf2(surf_info(index_mol,jj)%coord(2),&
                atom_info(m1,jj)%coord(3), thickness )

          !The following condition establish interface hydrogen bonds, 
          !which is the core of this method. 
          !Also check if the OH group is free OH
          IF ((condition1 .OR. condition2) .AND. H_info(k_H,jj)%IsFree .eqv. .True.) THEN
              nf(jj) = 1.0 
              freeoh_exist(jj) = .True.
              freeoh = freeoh + nf(jj) ! To calculate ave population of free OH over starting points for a OH group
          ENDIF
        END DO
        nfb(k) = freeoh 

        nfreeoh_exist(k) = nfreeoh
        tot_nfb = tot_nfb + freeoh
        tot_nfreeoh = tot_nfreeoh + nfreeoh_exist(k)

        !==============================================
        !Calcualte the correlation function C_freeOH(t)
        !==============================================
        if (nfb(k) > nfb_min) THEN
            DO mt = 0, nmo_effective-1 ! The time interval
                scalar = 0.0
                DO j = 1, nmo-nmo_effective, start_step ! How many steps? (nmo-nmo_effective-1)/start_step + 1 
                    tmp = nf(j) * nf(j+mt)
                    scalar = scalar + tmp 
                ENDDO
                corr_nf(mt+1) = corr_nf(mt+1) + scalar !sum_C_freeOH_k(t)
            ENDDO
        ENDIF
      ENDDO hLOOP
      nfb_per_frame = tot_nfb/REAL(nmo,rk)
      ave_nf = nfb_per_frame/REAL(nwat,rk) 
      !=========================================
      !Calculate the number of ever bonded pairs
      !=========================================
      n_freeoh = 0 
      DO k = 1, n_H
          IF (nf(k)>nfb_min) THEN
              n_freeoh = n_freeoh + 1      
          ENDIF
      ENDDO
      !==================================
      !Normalization of C_freeOH(t) step1
      !==================================
      corr_nf = corr_nf / (num_start_points * n_H)
      corr_nf = corr_nf / ave_nf

      !=====================================
      !Write the correlation
      !C_freeOH(t) for the iterfacial freeOH   
      !======================================
      char_thickness = nth(str(nint(thickness)),d_len)
      OPEN(10,file=trim(filename)//'_freeoh_acf_nf_' &
        //char_thickness//'.dat')
        DO i =1, nmo_effective
            WRITE(10,*) REAL(i-1, rk) * delta_t, corr_nf(i)
        ENDDO
        WRITE(6,*)'written in '//trim(filename)//&
                  '_freeoh_acf_nf_'//char_thickness//'.dat'
      close(10)
     !===========
     ! Print <nf>      
     !===========      
      OPEN(10,file=trim(filename)//'_freeoh_ave_nf_'//&
        char_thickness//'.dat')
        WRITE(10,*) 'Ave. No. free OHs:', nfb_per_frame
        WRITE(10,*) '<h>:', ave_nf
        WRITE(6,*)'written in '//trim(filename)//&
                  '_freeoh_ave_nf_'//char_thickness//'.dat'
      close(10)
      deallocate (nf,corr_nf,nfb)
      !==== to continue

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
        !ndx_3_list = h_ndx_list(ndx_1(k), ndx_2(k), pos_filename, nat, boxsize)
        ndx_3_list = (/ndx_H(2*idx_O1-1,2), ndx_H(2*idx_O1,2), ndx_H(2*idx_O2-1,2), ndx_H(2*idx_O2,2)/) 
        ! Calculate h(j)
        ! A LOOP on all frames
        TIME: DO jj = 1, nmo
          h(jj) = 0.0
          hb_exist(jj) = .False.

          ! Check if the pairs are located in one of the interfaces 
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
              HYDROGEN: DO j =1, 4
                  m3 = ndx_3_list(j)
                  r3 = (/atom_info(m3,jj)%coord(1),atom_info(m3,jj)%coord(2),&
                         atom_info(m3,jj)%coord(3) /)
                  r21 = dist2(r1, r2, boxsize) 
                  if (criterion == 1) THEN
                      r23 = dist2(r3, r2, boxsize) 
                      r13 = dist2(r3, r1, boxsize) 
                      pm = pmADH(r1,r2,r3,boxsize)  ! if H is bound to O2
                      pm_ = pmADH(r2,r1,r3,boxsize) ! if H is bound to O1
                      norm_rr = sqrt(r21*r23)
                      cosphi = pm/norm_rr
                      cosphi_= pm_/norm_rr
                      if ((r21 .lt. rooc) .and. ( (cosphi .gt. cosPhiC123) .or. &
                          (cosphi_  .gt. cosPhiC123) )                      &
                         ) THEN
                          h(jj) = 1.0 
                          hb_exist(jj) = .True.
                          qj = qj + h(jj) ! To calculate ave population of HB over all starting points for one pair of water                           
                          nqj = nqj + 1
                          EXIT ! if we know that two pair of molecule is bonded at frame jj, THEN we go to check the next frame (jj+1)
                      ENDIF
                  elseif (criterion == 2) THEN
                      !Follow the second cirterion of HB.
                      r31 = dist2(r1,r3,boxsize) 
                      r32 = dist2(r2,r3, boxsize) 
                      pm = pmAHD(r1,r2,r3,boxsize)
                      cosphi = pm/(sqrt(r31*r32))

                      !Follow the scond criterion of HB.
                      !-0.5 comes from cos(2pi/3)
                      if (r21 .lt. rooc .and. cosphi .lt. cosPhiC132) THEN 
                          h(jj) = 1.0 
                          hb_exist(jj) = .True.
                          qj = qj + h(jj) ! To calculate ave population of HB over all starting points for one pair of water                           
                          nqj = nqj + 1
                          EXIT ! If we know that two pair of molecule is bonded at frame jj, then we go to check the next frame (jj+1)
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
        !Calcualte the correlation function C_HB(t)
        !==========================================
        if (hb(k) > hb_min) THEN
            DO mt = 0, nmo_effective-1 ! The time interval
                scalar = 0.0
                DO j = 1, nmo-nmo_effective, start_step ! How many steps? (nmo-nmo_effective-1)/start_step + 1 
                    tmp = h(j)*h(j+mt)
                    scalar = scalar + tmp 
                ENDDO
                corr_h(mt+1) = corr_h(mt+1) + scalar !sum_C_k(t)
            ENDDO
        ENDIF
      ENDDO kLOOP   
      !deallocate(hb_exist,nhb_exist)
      hb_per_frame = tot_hb/REAL(nmo,rk)
      ave_h = hb_per_frame/REAL(nwat,rk) 
      !=========================================
      !Calculate the number of ever bonded pairs
      !=========================================
      n_bonded_pairs = 0 
      DO k = 1, nwat
          IF (hb(k)>hb_min) THEN
              n_bonded_pairs = n_bonded_pairs + 1      
          ENDIF
      ENDDO
      !==============================
      !Normalization of C_HB(t) step1
      !==============================
      corr_h = corr_h / (num_start_points*nwat)
      corr_h = corr_h / ave_h
      deallocate(x, y, z, ndx_1, ndx_2)          
      !===================================
      !Write the correlation
      !C_HB(t) for the iterfacial HB (ihb)    
      !===================================
      char_thickness = nth(str(nint(thickness)),d_len)
      OPEN(10,file=trim(filename)//'_wat_pair_hbacf_h_ihb_' &
        //char_thickness//'.dat')
        DO i =1, nmo_effective
            WRITE(10,*) REAL(i-1, rk) * delta_t, corr_h(i)
        ENDDO
        WRITE(6,*)'written in '//trim(filename)//&
                  '_wat_pair_hbacf_h_ihb_'//char_thickness//'.dat'
      close(10)
     !===========
     ! Print <h>      
     !===========      
      OPEN(10,file=trim(filename)//'_wat_pair_ave_h_ihb_'//&
        char_thickness//'.dat')
        WRITE(10,*) 'Ave. No. bonds:', hb_per_frame
        WRITE(10,*) '<h>:', ave_h
        WRITE(6,*)'written in '//trim(filename)//&
                  '_wat_pair_ave_h_ihb_'//char_thickness//'.dat'
      close(10)
      deallocate (h,corr_h,hb)
      deallocate (nf,corr_nf,nfb)
      END SUBROUTINE
