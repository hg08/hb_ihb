MODULE parameter_shared
!
! Purpose:
!   To declare data to share between routines.
IMPLICIT NONE
SAVE 
character(LEN=200) :: sampled_pos_filename
INTEGER, ALLOCATABLE, DIMENSION(:) :: sampled_movie
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: sampled_time, sampled_energy

! For griding
REAL(KIND=8), PARAMETER :: whish_size=0.5 ! Angstrom
INTEGER :: nb_divx, nb_divy, nb_divz, n_grid 
REAL(KIND=8) :: divx, divy, divz

! paramters of the interface
REAL(KIND=8) :: thickness ! the thickness of the instantaneous interfaces

END MODULE parameter_shared
