MODULE surf_traj
!
! Purpose: 
! To declear data related to the simulation and traj.
!
!========================================
! data         author         version
! ========     =======        ==========
! 2020-9-18    Gang Huang     Original
!========================================
!
IMPLICIT NONE

CONTAINS

  SUBROUTINE read_surf_traj(indx,nmo_start,nmo_end,ns,n_grid,n_samples)
    ! Purpose:
    ! To read info from the isosurface trajectory file (format: surf_traj_*.xyz)
    USE surf_module, ONLY: surf_info
    INTEGER :: i_sample, i_grid
    INTEGER, INTENT(IN) :: n_grid  ! KEEP
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns) KEEP
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER, INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    CHARACTER(LEN=4) :: head_char
    INTEGER :: y

    allocate(surf_info(n_grid,n_samples))
    i_sample = 1
    write(*,*) "read_surf_traj(): New total time steps (n_samples):", n_samples
    DO WHILE (i_sample < n_samples+1) ! +1 means i_sample can take the value of n_samples 
        read(indx, '(A4)') head_char
        PRE_CHECK:IF (head_char=="i = ") THEN
            BACKSPACE(UNIT=indx) ! Because I am not able to read other lines with the format '(A4,I8)', and have not find any good way, so I try to read it in '(A4)' first 
            read(indx, '(A4,I5)') head_char, y
            !read(indx, '(1X,A4,I5)') head_char, y
            CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-(nmo_start-1),ns) == 0) THEN
                !-------------------------------------------------------------------------------------------------------
                !NOTE: if use ' i = ', instead of 'i = ', it will be wrong!
                !IF (head_char==' i = ' .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-(nmo_start-1),ns) == 1) THEN
                !-------------------------------------------------------------------------------------------------------
                !WRITE(*,*)"read_traj():", head_char, y
                BACKSPACE(UNIT=indx) ! Because we have to read the whole line with ' i = ' line.
                READ(indx,*) ! skip one line in the unit=indx
                131 FORMAT (11X,2F13.6)
                inner: do i_grid= 1, n_grid
                  read (indx,131) surf_info(i_grid,i_sample)%coord(1), & 
                    surf_info(i_grid,i_sample)%coord(2)
                    !WRITE (*,131) surf_info(i_grid,i_sample)%coord(1), &
                    !surf_info(i_grid, i_sample)%coord(2)
                enddo inner
                i_sample = i_sample + 1 ! The position is important. It must be located before ENDIF MATCH
            ENDIF CHECK_HEAD
        ENDIF PRE_CHECK
    END DO

  END SUBROUTINE read_surf_traj

END MODULE surf_traj
