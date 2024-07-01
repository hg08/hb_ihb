MODULE water_module
! Purpose:
! To define the derived data type for water
USE atom_module

IMPLICIT NONE
TYPE :: water
  INTEGER :: molecule_id
  TYPE (atom):: oxygen
  TYPE (atom):: hydrogen1
  TYPE (atom):: hydrogen2
  REAL :: mass
END TYPE water

! The array water_info can be shared by subroutines  
TYPE(water), ALLOCATABLE, DIMENSION(:,:) :: water_info
END MODULE water_module
