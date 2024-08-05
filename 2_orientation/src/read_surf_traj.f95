SUBROUTINE read_surf_traj(surf_filename,nmo_start,nmo_end,ns,n_grid,n_samples,surf_info_fortran)
  ! Purpose:
  ! To read info from the isosurface trajectory file (format: surf_traj_*.xyz)
  ! In the v3 version, surf_info_fortran is an array of real values, instead of an array
  ! of a self-defined type.
  !=====================================
  ! data         author         version
  !=======     =======        ==========
  ! 2020-9-18    Gang Huang     Original
  !=====================================
  IMPLICIT NONE

  INTEGER :: ierror, i_sample, i_grid, ix, iy
  INTEGER,INTENT(IN) :: n_grid ! n_grid = nb_divx * nb_divy; KEEP
  INTEGER,INTENT(IN) :: n_samples ! n_samples = INT(nmo/ns); KEEP
  INTEGER,INTENT(IN) :: ns ! Get one sample from the trajectory every ns step.
  INTEGER,INTENT(IN) :: nmo_start, nmo_end ! To get the total number of moves
  REAL(kind=8),DIMENSION(2,n_samples,n_grid),INTENT(INOUT) :: surf_info_fortran
  CHARACTER(LEN=4) :: head_char, e
  character(LEN=200) :: surf_filename
  INTEGER :: y
  INTEGER :: skip_num ! The number of steps skipped before reading
  INTEGER, PARAMETER :: indx = 20

  !INITIALIZATION
  i_sample = 0
  skip_num = 0

  131 FORMAT (11X,2F13.6)
  open(indx,FILE=surf_filename,STATUS='OLD',ACTION='READ',IOSTAT=ierror)
  IF (nmo_start > 0) THEN
        DO i_sample = 0, nmo_start ! Skip
             read (indx,*)
             do i_grid = 1, n_grid
                read (indx,*)
             enddo
        END DO
  ENDIF

  DO i_sample = 0, n_samples - 1
        read(indx, *) head_char, e, y
        do i_grid = 1, n_grid
            read(indx, *, IOSTAT=ierror) ix, iy, surf_info_fortran(1,i_sample+1,i_grid), &
            surf_info_fortran(2,i_sample+1,i_grid)
        enddo 
  END DO   
  close(indx)

END SUBROUTINE read_surf_traj
