      SUBROUTINE ghbacf_interface_nf_pbc_format2(boxsize,delta_t0, &
          filename,pos_filename,n_samples,nat,ns,&
          criterion,atom_info,n_grid,divx,divy,divz,nb_divx,nb_divy,&
          nb_divz,thickness,surf_info)
      !========================================================================
      !Purpose: to obtain nf_HB(t): C_{freeOH}(t)= <nf(0)nf(t)>/<nf> for interfacial free OH
      !1)the correlation function is the average over OH groups of molecules. 
      !2)One can choose any time step to start, and end at any step for 
      !calculation, instead of the original time step.      
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
      USE atom_module, ONLY: atom
      USE water_molecule_types, ONLY: oxygen_atom, hydrogen_atom, O_info, H_info
      use module_ihb, ONLY: get_total_number_of_lines, &
                       hydrogen_ndx_list, &
                       distance2, &
                       pm_adh, pm_ahd, & 
                       grid_index, pair_in_surf1, pair_in_surf2,&
                       atom_in_surf1, atom_in_surf2, &
                       str, nth
      implicit none
      
      !==========
      !parameters
      !==========
      INTEGER, PARAMETER :: rk=8 ! local 
      INTEGER, PARAMETER :: d_len=1 ! for storing the length of the character which represents the thickness of the interface
      
      character(LEN=200), INTENT(INOUT) :: filename,pos_filename
      !character(LEN=200), INTENT(IN) :: list_filename
      INTEGER, INTENT(IN) :: criterion
      INTEGER, INTENT(IN) :: nat ! number of atoms
      INTEGER, INTENT(IN) :: ns
      INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
      REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: boxsize
      REAL(KIND=rk),INTENT(IN) :: delta_t0

      TYPE(atom), DIMENSION(nat,n_samples), INTENT(IN) :: atom_info
      
      REAL(KIND=rk), INTENT(IN) :: thickness ! the thickness of the instantaneous interfaces
      !REAL, PARAMETER :: whish_size=0.5 ! Angstrom
      INTEGER, INTENT(IN) :: nb_divx, nb_divy, nb_divz, n_grid 
      REAL(KIND=rk), INTENT(IN) :: divx, divy, divz
      REAL(KIND=rk),DIMENSION(2,n_grid,n_samples),INTENT(IN) :: surf_info

      !Local variables
      REAL(KIND=rk), PARAMETER :: rooc=12.25d0                 ! cutoff distance of rOO (3.5**2 )
      REAL(KIND=rk), PARAMETER :: cosPhiC123=0.866             !  phiC= 50 (degree).
      REAL(KIND=rk), PARAMETER :: cosPhiC132=-0.5            !  phiC132= 110 (degree).
      REAL(KIND=rk), PARAMETER :: cosPhiC123_freeOH=0.643             !  phiC= 50 (degree).
      REAL(KIND=rk), PARAMETER :: cosPhiC132_freeOH=-0.342            !  phiC132= 110 (degree).
      REAL(KIND=rk),PARAMETER :: h_min=0.5 ! condition for the existence of a h-bond for a step
      REAL(KIND=rk),PARAMETER :: max_time_for_corr = 12.0 ! Unit: ps.
      REAL(KIND=rk), PARAMETER :: hb_min=0.5 ! condition for the existence of h-bond for a pair of water molecules
      REAL(KIND=rk), PARAMETER :: nfb_min=0.5 ! condition for the existence of free OH for a  OH in water molecule
      REAL(KIND=rk) :: r13,cosphi,pm, cosphi_, pm_, norm_rr
      REAL(KIND=rk) :: r21,r31,r32,r23 ! For the second criterion of HB
      REAL(KIND=rk) :: qj,delta_t,hb_per_frame,ave_h
      REAL(KIND=rk), DIMENSION(3) :: r1, r2, r3 ! pbc 
      REAL(KIND=rk) :: freeoh, tot_nfb, nfb_per_frame, ave_nf
      INTEGER :: m1,m2,m3,mt,nqj,n_bonded_pairs
      INTEGER :: nfreeoh, tot_nfreeoh, n_freeoh
      INTEGER :: idx_O1
      INTEGER :: idx_H
      !REAL(KIND=rk), allocatable,DIMENSION (:) :: h,hb,corr_h
      REAL(KIND=rk),ALLOCATABLE,DIMENSION (:) :: nf, nfb, corr_nf
      REAL(KIND=rk), allocatable,DIMENSION (:,:) :: x,y,z
      INTEGER,ALLOCATABLE,DIMENSION(:) :: nfreeoh_exist
      INTEGER,ALLOCATABLE,DIMENSION(:) :: nf_exist
      !INTEGER, allocatable,DIMENSION(:) :: ndx_1,ndx_2,nhb_exist
      INTEGER, DIMENSION(4)   :: ndx_3_list
      REAL(KIND=rk)  :: scalar, sq, tmp 
      !LOGICAL,allocatable,DIMENSION (:)  :: hb_exist
      LOGICAL,ALLOCATABLE,DIMENSION (:) :: freeoh_exist
      INTEGER  :: nmo  ! nmo is not necessary, we set nmo = n_samples
      INTEGER :: nmo_effective, start_step, num_start_points
      INTEGER :: n_H ! number of OH groups; or num of H atoms
      INTEGER :: n_O ! number of O atoms
      INTEGER :: i,j,k,jj, bond 
      INTEGER :: k_H, k_O, k_O1, k_O2
      CHARACTER(LEN=d_len) :: char_thickness ! for saving the thickness in the files' names
      INTEGER :: index_mol
      INTEGER :: index_mol1, index_mol2
      LOGICAL :: condition1, condition2
      INTEGER :: order_type
      INTEGER :: tmp_index
      INTEGER :: num_limit ! Temperary varible, represents the largest number of an O atoms' neighbors.
      !==============
      !Initialization
      !==============
      scalar = 0.d0; sq = 0.d0;
      pm =0.d0; cosphi =0.d0
      r21 = 0.d0; r23 = 0.d0
      r31 = 0.d0; r13= 0.d0; r32 = 0.d0
      nfb_per_frame = 0.0; tot_nfb = 0.0
      hb_per_frame = 0.d0
      r1 = 0.d0; r2 = 0.d0; r3 = 0.d0
      nmo = n_samples
      ndx_3_list=0
      index_mol = 0
      index_mol1=0; index_mol2=0
      condition1=.FALSE.
      condition2=.FALSE.
      nmo_effective = 0
      start_step = 1
      order_type = 0

      n_H = nat * 2/3
      n_O = nat / 3
      allocate(nfb(n_H))
       
      delta_t=REAL(ns,rk)*delta_t0  ! unit: ps
      nmo_effective = nint(max_time_for_corr/delta_t) + 1
      start_step = nint((nmo_effective-1)/10.0) ! Start step of sliding window. Over-using rate is 1 - 1/5 = 4/5
      num_start_points = (nmo-nmo_effective-1)/start_step + 1
      ALLOCATE(x(nat,nmo))
      ALLOCATE(y(nat,nmo))
      ALLOCATE(z(nat,nmo))
      ALLOCATE(nf(nmo))
      allocate(nfreeoh_exist(n_H))
      !====================================
      ! Calculate <nf(0)nf(t)>/<nf>  
      ! Notice here <> is average over
      ! different OH groups in water molecules,
      ! and over all starting time points i
      ! with nf(i)=1.
      !====================================      
      !ALLOCATE(hb_exist(nmo))
      ALLOCATE(corr_nf(nmo))
      ALLOCATE(freeoh_exist(nmo))
      ! loop
      corr_nf(:)=0.d0 
      tot_nfb = 0
      tot_nfreeoh = 0
      nf(:) = 0
      nfb(:) = 0     
      nfreeoh_exist(:) = 0
      
      !! Check the O H arrange order.
      !if (atom_info(1,1)%atom_name == "O" .and. atom_info(2,1)%atom_name == "O") then
      !  ! O O ... H H H H ...
      !  order_type = 1
      !else
      !  ! O H H O H H ...
      !  order_type = 0
      !endif


      ! Calculate the correlation function
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
                condition1 = atom_in_surf1(surf_info(1,index_mol,jj),&
                      atom_info(m1,jj)%coord(3), thickness )
                !For surf 2 
                condition2 = atom_in_surf2(surf_info(2,index_mol,jj),&
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
      char_thickness = nth(str(thickness), d_len)
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

