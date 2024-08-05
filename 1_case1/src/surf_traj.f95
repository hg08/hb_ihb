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

  INTEGER :: ierror, i_sample, i_grid, ix, iy
  INTEGER,INTENT(IN) :: n_grid  ! KEEP
  INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns) KEEP
  INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
  INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
  REAL(KIND=8),DIMENSION(2,n_grid,n_samples),INTENT(INOUT) :: surf_info
  CHARACTER(LEN=4) :: head_char, e 
  CHARACTER(LEN=200) :: surf_filename
  INTEGER :: y
  INTEGER, PARAMETER :: indx = 20

  i_sample = 0 ! Frame number
  open(indx,FILE=surf_filename,STATUS='OLD',ACTION='READ',IOSTAT=ierror)
  IF (nmo_start > 0) THEN
        DO i_sample = 0, nmo_start - 1 ! Skip 
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
            read(indx, *, IOSTAT=ierror) ix, iy, surf_info(1,i_grid,i_sample+1), surf_info(2,i_grid,i_sample+1)
            !write(*,*) 'ix iy Surf1 Surf2', ix, iy, surf_info(1,i_grid,i_sample+1), surf_info(2,i_grid,i_sample+1)
        enddo 
  END DO 
  close(indx)
!write(*,*) "nmo_start, nmo_end, ns, n_grid, n_samples = ", nmo_start, nmo_end, ns, n_grid, n_samples
END SUBROUTINE read_surf_traj
