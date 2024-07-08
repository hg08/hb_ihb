SUBROUTINE read_surf_traj(surf_filename,nmo_start,nmo_end,ns,n_grid,n_samples,surf_info)
  ! Purpose:
  ! To read info from the isosurface trajectory file (format: surf_traj_*.xyz)
  ! In the v3 version, surf_info is an array of real values, instead of an array
  ! of a self-defined type.
  !=====================================
  ! data         author         version
  !=======     =======        ==========
  ! 2020-9-18    Gang Huang     Original
  !=====================================
  IMPLICIT NONE

  INTEGER :: ierror, i_sample, i_grid
  INTEGER,INTENT(IN) :: n_grid  ! KEEP
  INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns) KEEP
  INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
  INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
  REAL(KIND=8),DIMENSION(2,n_grid,n_samples),INTENT(INOUT) :: surf_info
  CHARACTER(LEN=4) :: head_char
  CHARACTER(LEN=200) :: surf_filename
  INTEGER :: y
  INTEGER, PARAMETER :: indx = 20


  OPEN(indx,FILE=surf_filename,STATUS='OLD',ACTION='READ',IOSTAT=ierror)
      ! TESTING
      !write(*,*) "read_surf_traj(): New total time steps (n_samples):", n_samples

      i_sample = 1
      !DO WHILE (i_sample < n_samples+1) ! +1 means i_sample can take the value of n_samples 
      DO WHILE (i_sample < n_samples) ! without '+1' means i_sample will not take the value of 'n_samples' 
          131 FORMAT (11X,2F13.6)
          READ(indx, '(A4)',IOSTAT=ierror) head_char
          PRE_CHECK:IF (head_char=="i = ") THEN
              BACKSPACE(UNIT=indx) ! Because I am not able to read other lines with the format '(A4,I8)', and have not find any good way, so I try to read it in '(A4)' first 
              READ(indx, '(A4,I5)',IOSTAT=ierror) head_char, y
              CHECK_HEAD: IF (head_char=="i = " .AND. (y>nmo_start .and. y<nmo_end+1) .AND. MOD(y-nmo_start,ns)==0 ) THEN
                  BACKSPACE(UNIT=indx) ! Because we have to read the whole line with ' i = ' line.
                  READ(indx,*) ! skip one line in the unit=indx
                  INNER: do i_grid= 1, n_grid
                    READ (indx,131,IOSTAT=ierror) surf_info(1,i_grid,i_sample), & 
                      surf_info(2,i_grid,i_sample)
                  ENDDO inner
                  i_sample = i_sample + 1 ! Located before ENDIF CHECK_HEAD
              ELSE
                  !WRITE(*,*) "None is satisfied."
                  i_sample = i_sample + 1 ! Located before ENDIF CHECK_HEAD 
              ENDIF CHECK_HEAD
          ENDIF PRE_CHECK
      END DO
  close(indx)

END SUBROUTINE read_surf_traj
