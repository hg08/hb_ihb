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
