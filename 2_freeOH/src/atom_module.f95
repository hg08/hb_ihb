MODULE atom_module
! Purpose:
!   To define the derived data type for atom
!
IMPLICIT NONE
TYPE :: atom
  CHARACTER(LEN=2) :: atom_name
  INTEGER :: atom_id 
  INTEGER :: host_id  ! For O atom in water, host_id = atom_id
  REAL(KIND=8) :: mass
  REAL(KIND=8),DIMENSION(3) :: coord 
END TYPE atom

! The array atom_info can be shared by subroutines  
TYPE(atom),ALLOCATABLE,DIMENSION(:,:) :: atom_info

END MODULE atom_module

module water_molecule_types
    implicit none

    ! Define the derived type for neighboring O atoms
    type :: oxygen_atom
        CHARACTER(LEN=2) :: atom_name
        INTEGER :: atom_id 
        INTEGER :: host_id  ! For O atom in water, host_id = atom_id
        REAL(KIND=8) :: mass
        REAL(KIND=8),DIMENSION(3) :: coord 
        integer, dimension(2) :: H_ids ! H1-O-H2 form a water molecule: H_ids(1) for H1, H_ids(2) for H2. 
        integer, dimension(12) :: self_indices_oxygen_neighbors 
        integer, dimension(12) :: indices_oxygen_neighbors 
        integer :: num_oxygen_neighbors
    end type oxygen_atom

    ! Define the derived type for H atoms
    type :: hydrogen_atom
        logical :: IsFree ! Whether the OH group is free
        CHARACTER(LEN=2) :: atom_name
        INTEGER :: atom_id 
        INTEGER :: host_id 
        REAL(KIND=8) :: mass
        REAL(KIND=8),DIMENSION(3) :: coord 
    end type hydrogen_atom

! The array O_info can be shared by subroutines  
TYPE(oxygen_atom), ALLOCATABLE, DIMENSION(:, :) :: O_info
! The array H_info can be shared by subroutines  
TYPE(hydrogen_atom), ALLOCATABLE, DIMENSION(:, :) :: H_info

end module water_molecule_types

