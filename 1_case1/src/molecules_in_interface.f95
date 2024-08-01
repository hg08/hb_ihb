      SUBROUTINE molecules_in_interface(n_samples,nat,arr1,arr2,atom_info,&
         n_grid,divx,divy,divz,nb_divx,nb_divy,nb_divz,thickness,&
         surf_info)
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
        USE module_ihb, ONLY: grid_index, mol_in_surf1, mol_in_surf2,&
                              str,nth,atom
        IMPLICIT NONE
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
        INTEGER :: jj,m,n1,n2 
        !To save the indices of the molecules for generating list file, we define a array for each time point (jj, in this code)
        INTEGER, DIMENSION(n_samples,nat), INTENT(INOUT) :: arr1, arr2
        TYPE(atom), DIMENSION(nat,n_samples), INTENT(IN) :: atom_info
        REAL(kind=8), DIMENSION(2,n_grid,n_samples), INTENT(IN) :: &
            surf_info

        !==============
        !Initialization
        !==============
        condition1=.FALSE.
        condition2=.FALSE.
        
        !initialize the arrays
        arr1 = 0
        arr2 = 0
       
        !=============
        !The main loop
        !=============      
        MOLECULE: DO jj =1, n_samples
          n1 = 0
          n2 = 0
          DO m=1,nat
            ! We use Oxygen atom to idenfy the water molecule
            IF (TRIM(atom_info(m,jj)%atom_name) == "O") THEN
              ! Check if the molecue is located in one of the interfaces 
              !write(*,*) "m", m
              !write(*,*) "jj", jj
              !write(*,*) "atom_info(m,jj)%coord(1)", atom_info(m,jj)%coord(1)
              !write(*,*) "atom_info(m,jj)%coord(2)", atom_info(m,jj)%coord(2)
              !write(*,*) "divx", divx
              !write(*,*) "divy", divy
              index_mol = grid_index(atom_info(m,jj)%coord(1), &
                  atom_info(m,jj)%coord(2),divx,divy,nb_divx,nb_divy) 
              !write(*,*) "index_mol", index_mol
              !For surf 1
              condition1 = mol_in_surf1(surf_info(1,index_mol,jj),&
                  atom_info(m,jj)%coord(3), thickness ) 
              !For surf 2 
              condition2 = mol_in_surf2(surf_info(2,index_mol,jj),&
                  atom_info(m,jj)%coord(3), thickness ) 

              !If the atom is in either surf1 or surf2
              IF (condition1) THEN
                  n1 = n1+1
                  !! Now write out the data into a 2-D array: arr1.
                  arr1(jj,n1)=m ! Jie: There should be two arrays, one for surf1, the other for surf2
              ELSE IF (condition2) THEN
                  n2 = n2+1
                  !! Now write out the data into a 2-D array: arr2.
                  arr2(jj,n2)=m 
              ENDIF
            ENDIF
          END DO

        ENDDO MOLECULE 
        
      END SUBROUTINE
