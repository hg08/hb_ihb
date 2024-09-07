      SUBROUTINE ghbacf_interface_nf_pbc_format2(boxsize, delta_t0, &
          filename, n_samples, nat, ns, &
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
      use tools, ONLY: get_nline => get_total_number_of_lines, &
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
      CHARACTER(LEN=200),INTENT(INOUT) :: filename
      !CHARACTER(LEN=200),INTENT(IN) :: list_filename_H
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
      REAL(KIND=rk) :: delta_t
      REAL(KIND=rk) :: freeoh, tot_nfb, nfb_per_frame, ave_nf 
      REAL(KIND=rk),DIMENSION(3) :: r1, r2, r3 ! pbc 
      INTEGER :: m1, m2, m3, mt, ns
      INTEGER :: nfreeoh, tot_nfreeoh, n_freeoh
      INTEGER :: idx_O1 ! The self index of O1 in all O atoms
      INTEGER :: idx_H ! The self Index of H in all H atoms
      REAL(KIND=rk),ALLOCATABLE,DIMENSION (:) :: nf, nfb, corr_nf
      REAL,ALLOCATABLE,DIMENSION (:,:) :: x, y, z
      INTEGER,ALLOCATABLE,DIMENSION(:) :: nfreeoh_exist
      INTEGER,ALLOCATABLE,DIMENSION(:) :: nf_exist
      !INTEGER,ALLOCATABLE,DIMENSION(:,:) :: ndx_H1, ndx_O1
      INTEGER,DIMENSION(4) :: ndx_3_list
      REAL(KIND=rk) :: scalar, tmp 
      LOGICAL,ALLOCATABLE,DIMENSION (:) :: freeoh_exist
      LOGICAL,ALLOCATABLE,DIMENSION (:) :: is_free
      INTEGER :: nmo ! nmo is not necessary, we set nmo = n_samples, because we DO not want to change too much
      INTEGER :: nmo_effective, start_step, num_start_points 
      INTEGER :: nwat ! number of water molecules
      INTEGER :: n_H ! number of OH groups; or num of H atoms
      INTEGER :: n_O ! number of O atoms
      INTEGER :: i, j, k, jj, bond 
      INTEGER :: k_H, k_O, k_O1, k_O2
      CHARACTER(LEN=d_len) :: char_thickness ! for saving the thickness in the files' names
      INTEGER :: index_mol1, index_mol2
      INTEGER :: index_mol
      LOGICAL :: condition1, condition2
      INTEGER :: tmp_index 
      INTEGER :: num_limit ! Temperary varible, represents the largest number of an O atoms' neighbors.
      !==============
      !Initialization
      !==============
      scalar = 0.0
      pm = 0.0; cosphi = 0.0
      r21 = 0.0; r23 = 0.0
      r31 = 0.0; r13 = 0.0; r32 = 0.0
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
      num_limit = 0

      !To obtain the total number of water pairs
      n_H = nat * 2/3 
      n_O = nat / 3
      allocate(nfb(n_H))          
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
      allocate(nf(nmo))
      allocate(nfreeoh_exist(n_H))
     !====================================
     ! Calculate <nf(0)nf(t)>/<h>  
     ! Notice here <> is not average over
     ! OH groups in water molecules,
     ! and over all starting time points i
     ! with nf(i)=1.
     !====================================      
      allocate(corr_nf(nmo))
      allocate(freeoh_exist(nmo))
      ! loop
      corr_nf(:) = 0.0      
      tot_nfb = 0.0
      tot_nfreeoh = 0
      nf(:) = 0.0 
      nfb(:) = 0.0
      nfreeoh_exist(:) = 0 

      !A loop over TOTAL Time Step. Purpose: Find Oxygen atoms' Oxygen neighbors.
      TTLOOP: DO jj = 1,nmo
      DO k_O1 = 1, n_O
          ! For a given O, we calculate the list of Oxygen neighbors of this O atom. 
          tmp_index = 0
          m1 = 3 * k_O1 - 2 ! Host O1, total index
          !write(*,*) "m1=", m1
          DO k_O2 = 1 , n_O
              if (k_O2 .ne. k_O1) then
                  m2 = 3 * k_O2 - 2 ! Total order of O2. For searching neighbors of Host O1
                  ! The total inices of H atoms related to O1 (Host O) and O2 (potential neighor O)
                  ! Check if O2 is a neighbor of O1
                  r1 = (/atom_info(m1,jj)%coord(1),atom_info(m1,jj)%coord(2),&
                     atom_info(m1,jj)%coord(3) /)
                  r2 = (/atom_info(m2,jj)%coord(1),atom_info(m2,jj)%coord(2),&
                         atom_info(m2,jj)%coord(3) /)
                  r21 = dist2(r1, r2, boxsize)
                  if ((r21 .lt. rooc)) THEN
                     !write(*,*) "m2=", m2
                     ! In the neighbor list of O1, add the total indices of O2.
                     O_info(k_O1,jj)%num_oxygen_neighbors = O_info(k_O1,jj)%num_oxygen_neighbors + 1 ! increase the number of O neighbors
                     tmp_index = O_info(k_O1,jj)%num_oxygen_neighbors
                     !write(*,*) "num_oxygen_neighbors for",k_O1,"-th O:", jj, "-th step", tmp_index
                     O_info(k_O1,jj)%self_indices_oxygen_neighbors(tmp_index) = k_O2
                     O_info(k_O1,jj)%indices_oxygen_neighbors(tmp_index) = m2 
                     O_info(k_O2, jj )%atom_id = m2 ! define total index of O2
                  endif
              endif          
          ENDDO
      ENDDO
      ENDDO TTLOOP

      !The purpose of this loop is to  define the free OH groups at each time step.
      tLOOP: DO jj = 1,nmo
          !===================================
          !The loop to find out free OH groups
          !===================================
          Ohost: DO k_O1 = 1, n_O ! self index of O1
              idx_O1 = k_O1 * 3 - 2 ! The total index of O1
              m1 = idx_O1 ! Indices of total index of the Host 
              OH: DO bond = 1,2
                  idx_H = O_info(k_O1,1)%H_ids(bond) ! 1: the first time step
                  m3 = idx_H 
                  if (bond == 1) then
                      k_H = (idx_H + 1) * 2/3 - 1
                  else
                      k_H = (idx_H ) * 2/3 
                  endif
                  if (m1 < 1) then
                      write(*,*) "Error: m1 must be positive!"
                  endif
                  if (O_info(k_O1,jj)%num_oxygen_neighbors .le. 0) then
                      H_info(k_H,jj)%IsFree = .True.
                  else
                      num_limit = O_info(k_O1,jj)%num_oxygen_neighbors
                      allocate(is_free(num_limit))
                      DO k_O2 = 1, num_limit 
                          m2 = O_info(k_O1,jj)%indices_oxygen_neighbors(k_O2) ! Indices of total index of the neighbor O2 
                          ! DONE: 一旦知此列表序列,对于任何一个氢原子,我们可以找到其宿主氧原子idx_O_host (H_info(:,:)%host_id)
                          ! For all tuples (idx_O1, idx_O_neighbor[1]), (idx_O1, idx_O_neighbor[2]), (idx_O1, idx_O_neighbor[3]),...
                          !Check whether the OH is free OH group respectively. Then based on these results, we can know that whether
                          !the OH is free OH or not. 
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
                              !norm_rr = sqrt(r21*r23) ! if H is bound to O2 (This is impossible here.)
                              norm_rr = sqrt(r21*r13) ! if H is bounded to O1
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
                      deallocate(is_free)
                  endif          
              ENDDO OH
          ENDDO Ohost
      ENDDO tLOOP

      OLOOP: DO k_O1 = 1, n_O ! Indices of self-indices of Host O
        freeoh = 0
        nfreeoh = 0 ! The number of freeoh-ed times for k-th OH group
        m1 = 3 * k_O1 -2 ! Total index of the Host (O1) 
        OHBOND: DO bond = 1, 2
            idx_H = O_info(k_O1,1)%H_ids(bond) ! 1: the fisrt time step
            if (bond == 1) then
                k_H = (idx_H + 1) * 2/3 - 1
            else
                k_H = (idx_H ) * 2/3
            endif
            time: DO jj = 1, nmo
                nf(jj) = 0.0
                freeoh_exist(jj) = .False.
                index_mol = grid_index(atom_info(m1,jj)%coord(1), &
                  atom_info(m1,jj)%coord(2),divx,divy,nb_divx)  

                ! Check if O1 is in one of the two interfaces
                !For surf 1
                condition1 = atom_in_surf1(surf_info(index_mol,jj)%coord(1),&
                      atom_info(m1,jj)%coord(3), thickness )
                !For surf 2 
                condition2 = atom_in_surf2(surf_info(index_mol,jj)%coord(2),&
                      atom_info(m1,jj)%coord(3), thickness )

                !The following condition establish interface free OH groups, 
                !ie., check if the OH group is free OH AND is in the interface.
                IF ((condition1 .OR. condition2) .AND. H_info(k_H,jj)%IsFree .eqv. .True.) THEN
                    nf(jj) = 1.0 
                    freeoh_exist(jj) = .True.
                    freeoh = freeoh + nf(jj) ! To calculate ave population of free OH over starting points for a OH group
                ENDIF

            ENDDO time
            nfb(k_H) = freeoh 
            nfreeoh_exist(k_H) = nfreeoh
            tot_nfb = tot_nfb + freeoh
            tot_nfreeoh = tot_nfreeoh + nfreeoh_exist(k_H)
            !==============================================
            !Calcualte the correlation function C_freeOH(t)
            !==============================================
            if (nfb(k_H) > nfb_min) THEN
                DO mt = 0, nmo_effective-1 ! The time interval
                    scalar = 0.0
                    DO j = 1, nmo-nmo_effective, start_step ! How many steps? (nmo-nmo_effective-1)/start_step + 1 
                        tmp = nf(j) * nf(j+mt)
                        scalar = scalar + tmp 
                    ENDDO
                    corr_nf(mt+1) = corr_nf(mt+1) + scalar !sum_C_freeOH_k(t)
                ENDDO
            endif
        ENDDO OHBOND
      ENDDO OLOOP
      nfb_per_frame = tot_nfb/REAL(nmo,rk)
      ave_nf = nfb_per_frame/REAL(n_H,rk) 
      !==========================================
      !Calculate the number of ever free OH groups
      !==========================================
      n_freeoh = 0 
      DO k = 1, n_H
          IF ( nf(k) > nfb_min ) THEN
              n_freeoh = n_freeoh + 1      
          ENDIF
      ENDDO
      write(*,*) "Number of ever free OH groups", n_freeoh
      !==================================
      !Normalization of C_freeOH(t) step1
      !==================================
      corr_nf = corr_nf / (num_start_points * n_H)
      corr_nf = corr_nf / ave_nf

      !=====================================
      !Write the correlation
      !C_freeOH(t) for the iterfacial freeOH   
      !======================================
      char_thickness = nth(str(nint(thickness)), d_len)
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
        WRITE(10,*) '<nf>:', ave_nf
        WRITE(6,*)'written in '//trim(filename)//&
                  '_freeoh_ave_nf_'//char_thickness//'.dat'
      close(10)
      deallocate (nf,corr_nf,nfb)
      END SUBROUTINE
