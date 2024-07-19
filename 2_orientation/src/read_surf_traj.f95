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

  INTEGER :: ierror, i_sample, i_grid
  INTEGER,INTENT(IN) :: n_grid ! n_grid = nb_divx * nb_divy; KEEP
  INTEGER,INTENT(IN) :: n_samples ! n_samples = INT(nmo/ns); KEEP
  INTEGER,INTENT(IN) :: ns ! Get one sample from the trajectory every ns step.
  INTEGER,INTENT(IN) :: nmo_start, nmo_end ! To get the total number of moves
  REAL(kind=8),DIMENSION(2,n_samples,n_grid),INTENT(INOUT) :: surf_info_fortran
  CHARACTER(LEN=4) :: head_char
  character(LEN=200) :: surf_filename
  INTEGER :: y
  INTEGER :: skip_num ! The number of steps skipped before reading
  INTEGER, PARAMETER :: indx = 20

  !INITIALIZATION
  i_sample = 1
  skip_num = 0

  131 FORMAT (11X,2F13.6)
  open(indx,FILE=surf_filename,STATUS='OLD',ACTION='READ',IOSTAT=ierror)
      ! TESTING
      write(*,*) "read_surf_traj(): New total time steps (n_samples):", n_samples
      !1.Determine How many lines have to be skiped
      skip_num = FLOOR(nmo_start/FLOAT(ns))
      !WRITE(*,*) "skip_num =", skip_num
      IF (skip_num > 0) THEN
          !SKIP skip_num lines
          DO i_sample = 1, skip_num
              read (indx,*) 
              inner0: do i_grid = 1, n_grid
                  read (indx,*) 
              enddo inner0
          END DO
          !READ AFTER SKIPPING LINES
          outer: DO i_sample = 1, n_samples ! i_sample can take the value of n_samples; WE WILL NOT USE DO WHILE LOOP ANY MORE IN THE FUTURE. 
              read(indx, '(A4,I5)',IOSTAT=ierror) head_char, y
              !WRITE(*, '(A24,A4,I5)') "PRE_CHECK head_char, y: ", head_char, y
              inner1: do i_grid = 1, n_grid
                  read (indx,131,IOSTAT=ierror) surf_info_fortran(1,i_sample,i_grid), & 
                      surf_info_fortran(2,i_sample,i_grid)
                  !WRITE (*,131) surf_info_fortran(1,i_sample,i_grid),surf_info_fortran(2,i_sample,i_grid)
              enddo inner1
          END DO outer
      ELSE
          !DIRECTLY READ
          outer2: DO i_sample = 1, n_samples ! i_sample can take the value of n_samples; WE WILL NOT USE DO WHILE LOOP ANY MORE IN THE FUTURE. 
              read(indx, '(A4,I5)',IOSTAT=ierror) head_char, y
              !WRITE(*, '(A29,A4,I5)') "PRE_CHEC_RIGHT head_char, y: ", head_char, y
              inner2b: do i_grid = 1, n_grid
                  read (indx,131,IOSTAT=ierror) surf_info_fortran(1,i_sample,i_grid), & 
                      surf_info_fortran(2,i_sample,i_grid)
                  !WRITE (*,131) surf_info_fortran(1,i_sample,i_grid),surf_info_fortran(2,i_sample,i_grid)
              enddo inner2b
          END DO outer2
      ENDIF
  close(indx)

END SUBROUTINE read_surf_traj
