      SUBROUTINE molecules_in_interface(n_samples,nat,arr,atom_info,&
         n_grid,divx,divy,divz,nb_divx,nb_divy,nb_divz,thickness,&
         surf_info_fortran)
        !=========================================================
        !Purpose: to obtain the indices of molecules in interfaces
        !         based on selection of molecules.
        !
        ! 2020-9-23
        ! Output: an 2D array--- indx_array 
        !         the first axis: time labels (jj)
        !         the 2nd axis: atom's labels (m) 
        !=========================================================
        !Modules
        use module_ihb, ONLY: grid_index, mol_in_surf1, mol_in_surf2,&
                              str,nth,atom
        implicit none
        !==========
        !parameters
        !==========
        INTEGER, PARAMETER :: rk= 8
        INTEGER, INTENT(IN) :: nat ! number of atoms
        INTEGER, INTENT(IN) :: n_samples ! steps of trajectory, ie, n_samples 
        INTEGER, INTENT(IN) :: nb_divx, nb_divy, nb_divz, n_grid 
        REAL(KIND=rk), INTENT(IN) :: divx, divy, divz
        REAL(KIND=rk), INTENT(IN) :: thickness
        !Local variables
        LOGICAL :: condition1,condition2  
        INTEGER :: index_mol
        INTEGER :: jj,m,n 
        !To save the indices of the molecules for generating list file, we define a array for each time point (jj, in this code)
        INTEGER, DIMENSION(n_samples,nat), INTENT(INOUT) :: arr
        TYPE(atom), DIMENSION(nat,n_samples), INTENT(IN) :: atom_info
        REAL(kind=8), DIMENSION(2,n_samples,n_grid), INTENT(IN) :: &
            surf_info_fortran

        !==============
        !Initialization
        !==============
        index_mol=0
        condition1=.FALSE.
        condition2=.FALSE.
        
        ! Use array instead of linked list, it may be faster. 
        !initialize the array
        arr = 0
       
        !=============
        !The main loop
        !=============      
        MOLECULE: DO jj =1, n_samples
          n = 0
          !WRITE(*,*) "NAT= ", nat, "[molecules_in_interface()]"
          !WRITE(*,*) "MOL1 nb_divx ", nb_divx
          DO m=1,nat
            ! We use Oxygen atom to idenfy the water molecule, THIS ASSUMPTION IS IMPORTANT
            IF (TRIM(atom_info(m,jj)%atom_name) == "O") THEN
              ! Check if the molecue is located in one of the interfaces 
              !WRITE(*,*) "[MOLecules_in_interface] nb_divx = ", nb_divx
              index_mol = grid_index(atom_info(m,jj)%coord(1), &
                  atom_info(m,jj)%coord(2),divx,divy,nb_divx, nb_divy) 
              !For surf 1
              !write(*,*)"TESTING :"
              !write(*,*) "thickness = ", thickness
              !write(*,*)"test b :"
              condition1 = mol_in_surf1(surf_info_fortran(1,jj,index_mol),&
                  atom_info(m,jj)%coord(3), thickness ) 
              !For surf 2 
              condition2 = mol_in_surf2(surf_info_fortran(2,jj,index_mol),&
                  atom_info(m,jj)%coord(3), thickness ) 
              !WRITE(*,*) condition1, condition2 

              !If the atom is in either surf1 or surf2
              !IF (condition1 .OR. condition2) THEN
              IF (condition2) THEN
                  n = n+1
                  !! Now write out the data into a 2-D array: arr.
                  arr(jj,n)=m
                  !WRITE(*,*) "MOL Index of atom:", m
                  !WRITE(*,*) "MOL Atom name:", atom_info(m,jj)%atom_name
              ENDIF
            ENDIF
          END DO

        ENDDO MOLECULE 
        
      END SUBROUTINE 
