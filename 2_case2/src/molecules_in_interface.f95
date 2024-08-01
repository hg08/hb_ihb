      SUBROUTINE molecules_in_interface(nmo, nat, arr)
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
      use tools, ONLY: get_nwat => get_total_numer_of_lines, &
                       dist2 => distance2, &
                       pmADH, pmAHD, & 
                       grid_index, mol_in_surf1, mol_in_surf2,&
                       str, nth
      USE atom_module 
      USE parameter_shared
      USE traj_format2 
      USE surf_traj ! For reading data from isosurface trajectory
      USE surf_module, ONLY: surf_info
 
      implicit none
      
      !==========
      !parameters
      !==========
      INTEGER, INTENT(IN) :: nat ! number of atoms
      INTEGER, INTENT(IN) :: nmo ! steps of trajectory
      integer,parameter :: rk=4 ! local 

      !Local variables
      INTEGER :: jj,m,n 
      INTEGER :: index_mol
      LOGICAL :: condition1,condition2  

      !To save the indices of the molecules for generating list file, we define a array for each time point (jj, in this code)
      INTEGER, DIMENSION(nmo,nat), INTENT(INOUT) :: arr
      
      !==============
      !Initialization
      !==============
      index_mol=0
      condition1=.FALSE.
      condition2=.FALSE.
      
      ! Use array instead of linked list, it may be faster. 
      !allocate(arr(nmo,nat))

      !initialize the array
      arr = 0

      !=============
      !The main loop
      !=============      
      MOLECULE: do jj =1, nmo
        n = 0
        DO m = 1, nat
          ! We use Oxygen atom to idenfy the water molecule
          IF (TRIM(atom_info(m,jj)%atom_name) == "O") THEN
            ! Check if the molecue is located in one of the interfaces 
            index_mol = grid_index(atom_info(m,jj)%coord(1), &
                atom_info(m,jj)%coord(2),divx,divy,nb_divy) 
            !WRITE(*,*) "index_mol = ",index_mol

            !For surf 1
            condition1 = mol_in_surf1(surf_info(index_mol,jj)%coord(1),&
                atom_info(m,jj)%coord(3), thickness ) 

            !For surf 2 
            condition2 = mol_in_surf2(surf_info(index_mol,jj)%coord(2),&
                atom_info(m,jj)%coord(3), thickness ) 
            !WRITE(*,*) condition1, condition2 

            !If the atom is in either surf1 or surf2
            IF (condition1 .OR. condition2) THEN
                n = n+1
                !! Now write out the data into a -D array: arr.
                arr(jj,n)=m
                !WRITE(*,*) "Index of atom:", m
                !WRITE(*,*) "Atom name:", atom_info(m,jj)%atom_name
            ENDIF
          ENDIF
        END DO

      enddo MOLECULE 
      
      !!TEST: print parts of the array
      !WRITE (*,*) "n:",n
      !DO m = 1,nat
      !  WRITE (*,*) "index array: step 1", arr(3,m)
      !  WRITE (*,*) "index array: step 2", arr(4,m)
      !ENDDO 
      !!TEST: print the whole array
      !WRITE (*,*) "index array:", arr

      END SUBROUTINE 
