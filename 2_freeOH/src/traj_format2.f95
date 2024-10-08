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

    !write(*,*) 'In function sampling_number: nmo_end = ', nmo_end
    ! no. of samples = INT({no. of moves}/ns)
    positive: IF (nmo_end <0 .OR. nmo_start < 0 .OR. ns <0) THEN
      write(*,*) 'Please enter non-negative values for the ns, starting step and ending step.'
    ELSE IF (nmo_end < nmo_start) THEN
      write(*,*) 'Please note that starting step shoud not larger than  ending step.'
      ELSE IF (ns == 0) THEN
        sampling_number = nmo_end-(nmo_start-1)
      ELSE
        sampling_number = FLOOR(FLOAT(nmo_end-(nmo_start-1))/FLOAT(ns))
      END IF positive
!    ELSE IF (ns ==0) THEN
!      sampling_number = nmo_end-(nmo_start-1)
!    ELSE IF (nmo_end-(nmo_start-1) <= ns) THEN
!      sampling_number = INT((nmo_end-(nmo_start-1))/ns + 1)
!    ELSE IF (nmo_end-(nmo_start-1) > ns) THEN
!      sampling_number = INT((nmo_end-(nmo_start-1))/ns)
!    END IF positive
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
      !write(*,*) 'the step:', imovie
      inner: do iatom= 1,nat
        read (indx,*) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
          atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
        !WRITE (*,*) & 
        !atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), &
        !atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
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
    USE parameter_shared, ONLY: sampled_movie, sampled_time
    INTEGER :: iatom, i_sample
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
        PRE_CHECK:IF (head_char=="i = ") THEN
            BACKSPACE(UNIT=indx) ! Because I am not able to read other lines with the format '(A4,I8)', and have not find any good way, so I try to read it in '(A4)' first 
            !read(indx, '(A4,I8)') head_char, y
            read(indx, '(1X,A4,I8)') head_char, y
            CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-nmo_start,ns) == 0) THEN
            !Jie: For special case of ns=1, MOD(y-(nmo_start-1),ns) is always 0. Hence, it needs to be checked separately. 
            !Use ‘&’ to continue the line to avoid Fortran maximum line length of 132 characters. 
            !CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. &
            !                (ns == 1 .or. MOD(y-(nmo_start-1),ns) == 1))  THEN 
                BACKSPACE(UNIT=indx) ! Because we have to read the whole line with ' i = ' line.
                read(indx,130) sampled_movie(i_sample), sampled_time(i_sample)
                130 FORMAT (5X,I8,9X,F12.3,6X,F20.10)
                !131 FORMAT (A4,3F20.10)
                inner: do iatom= 1,nat
                  read (indx,*) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
                    atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
                  if (atom_info(iatom, i_sample)%atom_name == "O") THEN
                      atom_info(iatom, i_sample)%mass = 16.00
                  elseif (atom_info(iatom, i_sample)%atom_name == "H") THEN
                      atom_info(iatom, i_sample)%mass = 1.00
                  endif
                enddo inner
                i_sample = i_sample + 1 !The position is important. It must be located before ENDIF MATCH
            ENDIF CHECK_HEAD
        ENDIF PRE_CHECK
    END DO
  END SUBROUTINE read_traj

  SUBROUTINE read_traj_O_H(indx,nmo_start,nmo_end,ns,nat,n_samples)
    ! Purpose:
    ! To read info from the trajectory file (format: ***.xyz)
    ! The data are saved in atom_info, O_info and H_info.
    ! to READ data starting from a pattern-matched line.
    USE atom_module, ONLY: atom_info
    USE water_molecule_types, ONLY: O_info, H_info
    USE parameter_shared, ONLY: sampled_movie, sampled_time
    INTEGER :: iatom,i_sample, i_O, i_H
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER, INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    CHARACTER(LEN=4) :: head_char
    INTEGER :: y

    allocate(atom_info(nat,n_samples))
    allocate(O_info(nat/3,n_samples))
    allocate(H_info(2*nat/3,n_samples))
    i_sample = 1
    write(*,*) "read_traj(): New total time steps (n_samples):", n_samples
    DO WHILE (i_sample < n_samples+1) ! +1 means i_sample can take the value of n_samples 
        !read(indx, '(A4)') head_char
        read(indx, '(1X,A4)') head_char  ! for some other format, one can use this format
        PRE_CHECK:IF (head_char=="i = ") THEN
            BACKSPACE(UNIT=indx)
            read(indx, '(1X,A4,I8)') head_char, y
            CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-nmo_start,ns) == 0) THEN
            !Jie: For special case of ns=1, MOD(y-(nmo_start-1),ns) is always 0. Hence, it needs to be checked separately. 
                BACKSPACE(UNIT=indx) ! Because we have to read the whole line with ' i = ' line.
                read(indx,130) sampled_movie(i_sample), sampled_time(i_sample)
                130 FORMAT (5X,I8,9X,F12.3,6X,F20.10)
                !131 FORMAT (A4,3F20.10)
                i_O = 0
                i_H = 0
                inner: do iatom= 1,nat
                  read (indx,*) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
                    atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
                  if (atom_info(iatom, i_sample)%atom_name == "O") THEN
                      i_O = i_O + 1
                      atom_info(iatom, i_sample)%mass = 16.00
                      atom_info(iatom, i_O)%atom_name = "O"
                      O_info(i_O,i_sample)%coord = atom_info(iatom,i_sample)%coord
                      O_info(i_O,i_sample)%atom_id = iatom
                      O_info(i_O,i_sample)%host_id = iatom
                      O_info(i_O,i_sample)%H_ids(1) = iatom + 1
                      O_info(i_O,i_sample)%H_ids(2) = iatom + 2
                      O_info(i_O,i_sample)%num_oxygen_neighbors = 0
                      O_info(i_O,i_sample)%self_indices_oxygen_neighbors = 0
                      O_info(i_O,i_sample)%indices_oxygen_neighbors = 0
                  elseif (atom_info(iatom, i_sample)%atom_name == "H") THEN
                      i_H = i_H + 1
                      atom_info(iatom, i_sample)%mass = 1.00
                      atom_info(i_H, i_sample)%atom_name = "H"
                      H_info(i_H, i_sample)%coord = atom_info(iatom,i_sample)%coord
                      H_info(i_H, i_sample)%atom_id = iatom
                      ! Assgining host_id. 3: water molecule has 3 atoms.
                      if ( mod(iatom,3) == 2 ) then
                          H_info(i_H, i_sample)%host_id = iatom - 1 
                      else if ( mod(iatom,3) == 0 ) then
                          H_info(i_H, i_sample)%host_id = iatom - 2 
                      else
                         write(*,*) "The order of O and H is not `O H H O H H...`."
                      endif
                      H_info(i_H, i_sample)%IsFree = .False. 
                  endif
                enddo inner
                i_sample = i_sample + 1 !The position is important. It must be located before ENDIF MATCH
            ENDIF CHECK_HEAD
        ENDIF PRE_CHECK
    END DO
  END SUBROUTINE read_traj_O_H

  SUBROUTINE skip_lines(indx, i_input)
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
