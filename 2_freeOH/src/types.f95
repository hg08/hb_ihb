  MODULE types
  ! Derived data type to stroe array
  TYPE :: real_array
    CHARACTER(2) :: atom_name
    INTEGER :: time_step
    REAL(8), DIMENSION(3) :: coord_array
    TYPE (real_array), POINTER :: p
  END TYPE real_array
  ! Derived data type to store water molecules around a center atom
  TYPE :: water_array
    INTEGER :: time_step
    CHARACTER(2) :: center_atom_name
    REAL(8), DIMENSION(9) :: coord_water_array
    REAL(8) :: d_RO
    REAL(8), DIMENSION(2) :: theta 
    TYPE (water_array), POINTER :: p
  END TYPE water_array
  END MODULE types
