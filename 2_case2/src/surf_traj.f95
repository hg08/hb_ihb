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
    INTEGER :: ierror, i_sample, i_grid, ix, iy
    INTEGER, INTENT(IN) :: n_grid  ! KEEP
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns) KEEP
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER, INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    CHARACTER(LEN=4) :: head_char,e
    INTEGER :: y

    ALLOCATE(surf_info(n_grid,n_samples))
    i_sample = 0
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
          !write(*,*), head_char, e, y
          do i_grid = 1, n_grid
              !read(indx, *, IOSTAT=ierror) ix, iy, surf_info(1,i_grid,i_sample+1), surf_info(2,i_grid,i_sample+1)
              read (indx,*) ix, iy, surf_info(i_grid,i_sample)%coord(1), surf_info(i_grid,i_sample)%coord(2)
              write (*,*) surf_info(i_grid,i_sample)%coord(1), surf_info(i_grid,i_sample)%coord(2)
          enddo 
    END DO 
    close(indx)
  END SUBROUTINE read_surf_traj

END MODULE surf_traj
