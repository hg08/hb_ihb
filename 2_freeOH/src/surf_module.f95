MODULE surf_module
! Purpose:
!   To define the derived data type for isosurfaces of an interfacial system
IMPLICIT NONE
TYPE :: surf
  INTEGER :: grid_index_x, grid_index_y 
  REAL(kind=8), DIMENSION(2) :: coord ! upper and lower location of the isosurface two real values.
END TYPE surf

! The array surf_info can be shared by subroutines  
TYPE(surf), ALLOCATABLE, DIMENSION(:,:) :: surf_info
END MODULE surf_module
