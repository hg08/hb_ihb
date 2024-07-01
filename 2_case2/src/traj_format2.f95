MODULE traj_format2
!
! Purpose: 
! To declear data related to the simulation and traj.
IMPLICIT NONE
CONTAINS
  INTEGER FUNCTION sampling_number(nmo_start,nmo_end,ns)
    !
    ! Purpose:
    !  To calculate the total numbers of samples one want to include in their analysis.
    ! Data dictionary
    INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves

    write(*,*) 'In function sampling_number: nmo_end = ', nmo_end
    ! no. of samples = INT({no. of moves}/ns)
    positive: IF (nmo_end <0 .OR. nmo_start < 0 .OR. ns <0) THEN
      write(*,*) 'Please enter non-negative values for the ns, starting step and ending step.'
    ELSE IF (nmo_end < nmo_start) THEN
      write(*,*) 'Please note that starting step shoud not larger than  ending step.'
    ELSE IF (ns ==0) THEN
      sampling_number = nmo_end-(nmo_start-1)
    ELSE IF (nmo_end-(nmo_start-1) <= ns) THEN
      sampling_number = INT((nmo_end-(nmo_start-1))/ns + 1)
    ELSE IF (nmo_end-(nmo_start-1) > ns) THEN
      sampling_number = INT((nmo_end-(nmo_start-1))/ns)
    END IF positive
  END FUNCTION sampling_number

  SUBROUTINE read_traj_v1(indx,nmo_start,nmo_end,ns,nat,n_samples)
    ! Purpose:
    ! To read info from the trajectory file (format: ***.xyz)
    ! to READ data starting from a pattern-matched line.
    USE atom_module
    USE parameter_shared

    INTEGER :: iatom, imovie, i_sample
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    !TYPE(atom), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: atom_info
    allocate(atom_info(nat,n_samples))
    i_sample = 1
    outer: DO imovie=1,nmo_end-nmo_start
      read(indx,*)!Neglect data of this line
      read(indx,120) sampled_movie(i_sample), sampled_time(i_sample), sampled_energy(i_sample)
      120 FORMAT (5X,I8,9X,F12.3,6X,F20.10)
      write(*,*) 'the step:', imovie
      inner: do iatom= 1,nat
        read (indx,*) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
          atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
        WRITE (*,*) & 
        atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), &
        atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
      enddo inner
      ! ns should be non-negative
      IF (ns < 0) THEN
        WRITE(*,*) 'Please note that ns should be non-negative.'
      ELSE IF (ns == 0) THEN
        CYCLE
      ELSE
        CALL skip_lines(indx, (nat+2)*(ns-1)) ! use CALL to run a SUBROUTINE
      ENDIF 

      i_sample = i_sample + 1

      ! To check if the sampling is finished
      check: IF (i_sample > n_samples) THEN 
        WRITE (*,*)'The total number of sample points are: ', n_samples
        EXIT
      END IF check
    ENDDO outer

  END SUBROUTINE read_traj_v1

  SUBROUTINE read_traj(indx,nmo_start,nmo_end,ns,nat,n_samples)
    ! Purpose:
    ! To read info from the trajectory file (format: ***.xyz)
    ! to READ data starting from a pattern-matched line.
    USE atom_module, ONLY: atom_info
    USE parameter_shared, ONLY: sampled_movie, sampled_time, sampled_energy
    INTEGER :: iatom,i_sample
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER, INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    CHARACTER(LEN=4) :: head_char
    INTEGER :: y

    allocate(atom_info(nat,n_samples))
    i_sample = 1
    write(*,*) "read_traj(): New total time steps (n_samples):", n_samples
    DO WHILE (i_sample < n_samples+1) ! +1 means i_sample can take the value of n_samples 
        !read(indx, '(A4)') head_char
        read(indx, '(1X,A4)') head_char  ! for some other format, one can use this format
        !------------------------------------------ 
        !TEST
        !WRITE(*,*)"head_char AND y:", head_char, y
        !------------------------------------------ 
        PRE_CHECK:IF (head_char=="i = ") THEN
            !WRITE(*,*)"Debug: PRE_CHECK head_char:", head_char
            !WRITE(*,*)"Debug: PRE_CHECK i_sample:",  i_sample
            BACKSPACE(UNIT=indx) ! Because I am not able to read other lines with the format '(A4,I8)', and have not find any good way, so I try to read it in '(A4)' first 
            !read(indx, '(A4,I8)') head_char, y
            read(indx, '(1X,A4,I8)') head_char, y
            !WRITE(*,*)"Debug: PRE_CHECK y:", y 
            !WRITE(*,*)"Debug: PRE_CHECK nmo_start-1:", nmo_start-1 
            !WRITE(*,*)"Debug: PRE_CHECK nmo_end+1:", nmo_end+1 
            !WRITE(*,*)"Debug: PRE_CHECK MOD(y-(nmo_start-1),ns):", MOD(y-(nmo_start-1),ns)
            WRITE(*,*)"Debug: PRE_CHECK ns:", ns
            !CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-(nmo_start-1),ns) == 1) THEN
            !Jie: For special case of ns=1, MOD(y-(nmo_start-1),ns) is always 0. Hence, it needs to be checked separately. 
            !Use ‘&’ to continue the line to avoid Fortran maximum line length of 132 characters. (Stupid Fortran!)
            CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. &
                            (ns == 1 .or. MOD(y-(nmo_start-1),ns) == 1))  THEN 
                !-------------------------------------------------------------------------------------------------------
                !NOTE: if use ' i = ', instead of 'i = ', it will be wrong!
                !IF (head_char==' i = ' .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-(nmo_start-1),ns) == 1) THEN
                !-------------------------------------------------------------------------------------------------------
                WRITE(*,*)"read_traj():", head_char, y
                BACKSPACE(UNIT=indx) ! Because we have to read the whole line with ' i = ' line.
                read(indx,130) sampled_movie(i_sample), sampled_time(i_sample), sampled_energy(i_sample)
                130 FORMAT (5X,I8,9X,F12.3,6X,F20.10)
                !130 FORMAT (4X,I8,9X,F12.3,6X,F20.10)
                131 FORMAT (A4,3F20.10)
                inner: do iatom= 1,nat
                  !read (indx,131) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
                    !atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
                  read (indx,*) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
                    atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
                  if (atom_info(iatom, i_sample)%atom_name == "O") THEN
                      atom_info(iatom, i_sample)%mass = 16.00
                  elseif (atom_info(iatom, i_sample)%atom_name == "H") THEN
                      atom_info(iatom, i_sample)%mass = 1.00
                  elseif (atom_info(iatom, i_sample)%atom_name == "N") THEN
                      atom_info(iatom, i_sample)%mass = 14.00 
                  elseif (atom_info(iatom, i_sample)%atom_name == "Li") THEN
                      atom_info(iatom, i_sample)%mass = 6.94
                  elseif (atom_info(iatom, i_sample)%atom_name == "Na") THEN
                      atom_info(iatom, i_sample)%mass = 22.99
                  elseif (atom_info(iatom, i_sample)%atom_name == "K") THEN
                      atom_info(iatom, i_sample)%mass = 39.10
                  endif
                  !WRITE (*,131) & 
                  !atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), &
                  !atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
                enddo inner
                i_sample = i_sample + 1 !The position is important. It must be located before ENDIF MATCH
            ENDIF CHECK_HEAD
        ENDIF PRE_CHECK
    END DO
      ! To check if the sampling is finished
      !check: IF (i_sample > n_samples) THEN 
      !  WRITE (*,*)'The total number of sample points are: ', n_samples
      !END IF check
  END SUBROUTINE read_traj

  SUBROUTINE skip_lines(indx, i_input)
    !
    ! Purpose: 
    ! To skip lines when read data from the input
    IMPLICIT NONE
    INTEGER :: i
    INTEGER,INTENT(IN) :: i_input,indx
    do i=1,i_input
       read(indx,*) !Neglect (nat+2)*(ns-1) lines
    enddo    
  END SUBROUTINE skip_lines

END MODULE traj_format2
