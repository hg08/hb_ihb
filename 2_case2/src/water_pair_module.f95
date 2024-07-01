MODULE water_pair_module
!
! Purpose:
! To define the derived data type for a water pair

USE water_module
IMPLICIT NONE
TYPE :: water_pair
  INTEGER :: water_pair_id
  INTEGER :: head_id ! For example, 2
  INTEGER :: tail_id ! For example, 1
  TYPE (water) :: head ! water molecule
  TYPE (water) :: tail  !water molecule
  LOGICAL :: h_state 
  REAL :: h ! Hydrogen bond population operator 1 or 0 
  REAL :: hb_length 
END TYPE water_pair

!The array ```water_pair_info``` can be shared by subroutines  
TYPE(water_pair), ALLOCATABLE, DIMENSION(:,:) :: water_pair_info
END MODULE water_pair_module  
