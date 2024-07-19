SUBROUTINE sample_and_recenter_format2(pos_filename,nmo_start,nmo_end,nat,ns,n_samples,boxsize, &
           sampled_pos_filename,sampled_movie,sampled_time,sampled_energy, &
           nb_divx,nb_divy,nb_divz,n_grid,divx,divy,divz,whish_size,atom_info)
  !2018/12/27
  ! 0) Purpose:     
  ! 0a) The subroutine sample_and_recenter.f95 reduces the size of the trajectory and recenter each step of the trajectory. 
  ! 0b) Date     Programmer     version
  ! 20.d0-9-16   Gang Huang     1
  ! Declaration
  USE module_ihb, ONLY: atom,read_traj
      
  IMPLICIT NONE
  !==========
  !Parameters
  !==========
  integer, parameter :: rk=8

  ! Varables
  character(LEN=200), INTENT(IN) :: pos_filename
  INTEGER, INTENT(IN) :: nmo_start  
  INTEGER, INTENT(IN) :: nmo_end  
  INTEGER, INTENT(IN) :: nat ! number of atoms or nb_atoms
  INTEGER, INTENT(IN) :: ns ! Get one sample from the trajectory every ns step.
  INTEGER, INTENT(IN) :: n_samples !n_samples = INT(nmo/ns)
  real(kind=rk), dimension(3), INTENT(INOUT) :: boxsize
  REAL(KIND=rk), INTENT(IN) :: whish_size ! Angstrom
  character(LEN=200), INTENT(INOUT):: sampled_pos_filename
  INTEGER, DIMENSION(n_samples),INTENT(INOUT) :: sampled_movie
  REAL(KIND=rk), DIMENSION(n_samples), INTENT(INOUT) :: sampled_time, sampled_energy
  INTEGER, INTENT(INOUT) :: nb_divx, nb_divy, nb_divz, n_grid 
  REAL(KIND=rk), INTENT(INOUT) :: divx, divy, divz
  TYPE(atom), DIMENSION(nat,n_samples),INTENT(INOUT) :: atom_info
  ! Local varables
  INTEGER :: iatom,imovie,i
  INTEGER :: num_wat_pairs 
  INTEGER, PARAMETER :: n_axis = 2  ! n_axis=2 means z-axis is the normal axis.
  REAL(KIND=rk) :: sum_mass
  REAL(KIND=rk) :: center_pos 
  
  ! Initialization
  iatom = 0; imovie = 0; i = 0; num_wat_pairs = 0
  nb_divx = nint(boxsize(1)/whish_size) ! round the argument to the nearest integer.
  !WRITE(*,*) "SAMPLE nb_divx = ", nb_divx 
  !WRITE(*,*) "SAMPLE whish_size = ", whish_size
  nb_divy = nint(boxsize(2)/whish_size) ! round the argument to the nearest integer.
  nb_divz = nint(boxsize(3)/whish_size) ! round the argument to the nearest integer.
  divx = boxsize(1)/REAL(nb_divx,rk)
  divy = boxsize(2)/REAL(nb_divy,rk)
  divz = boxsize(3)/REAL(nb_divz,rk)
  n_grid = nb_divx * nb_divy
  !WRITE(*,*)"n_grid=", n_grid
  !WRITE(*,*)"SAMPLE2 nb_divx = ", nb_divx
  !For pure water system:
  !(TODO)Formula may have to added for the cases of other systems
  num_wat_pairs = (nat/3)*(nat/3-1)/2

  !=======================
  !read in trajectory file 
  !=======================
  open(10,file=trim(pos_filename))
  ! Now starting read data
  !write(*,*) "Total number of atoms: ", nat
  CALL read_traj(10,nmo_start,nmo_end,ns,nat,n_samples,sampled_movie,sampled_time,sampled_energy,atom_info) 
  close(10)
  write(6,*) 'End of trajectory reading.'

  !=============
  !write in file
  !=============
  sampled_pos_filename = 'recentered_traj_pos_sampled.xyz'
  open(10,file=sampled_pos_filename)

  step: DO i=1,n_samples
    sum_mass = 0.d0
    center_pos = 0.d0

    write (10,'(I8)') nat
    WRITE(10,100) ' i = ',i-1,', time = ',sampled_time(i),', E = ',sampled_energy(i)
    100 FORMAT (A5,I8,A9,F12.3,A6,F20.10)
    !130 FORMAT (5X,I8,9X,F12.3,6X,F20.10)
  
    DO iatom=1,nat
        if (n_axis==0) THEN
            center_pos=center_pos+atom_info(iatom,i)%coord(1)*atom_info(iatom,i)%mass
        elseif (n_axis==1) THEN 
            center_pos=center_pos+atom_info(iatom,i)%coord(2)*atom_info(iatom,i)%mass
        else
            center_pos=center_pos+atom_info(iatom,i)%coord(3)*atom_info(iatom,i)%mass
        endif
        sum_mass=sum_mass+atom_info(iatom,i)%mass
    ENDDO
    center_pos=center_pos/sum_mass
    
    ! Recenter (along the normal axis only) the atoms according to the mass center
    DO iatom=1,nat
        IF (n_axis == 0) THEN
            atom_info(iatom,i)%coord(1) = atom_info(iatom,i)%coord(1) - center_pos + boxsize(1)/2
        ELSEIF (n_axis == 1) THEN
            atom_info(iatom,i)%coord(2) = atom_info(iatom,i)%coord(2) - center_pos + boxsize(2)/2
        ELSE 
            atom_info(iatom,i)%coord(3) = atom_info(iatom,i)%coord(3) - center_pos + boxsize(3)/2
        ENDIF
    ENDDO     

    ! Format for writing in recentered 
    !200 FORMAT (1X,A3,3F20.10)
    131 FORMAT (A4,3F20.10)
     
    DO iatom = 1, nat

        atom_info(iatom,i)%coord(1)=mod(atom_info(iatom,i)%coord(1),boxsize(1))
        if (atom_info(iatom,i)%coord(1) < 0) THEN
           atom_info(iatom,i)%coord(1) = atom_info(iatom,i)%coord(1) + boxsize(1) 
        endif
        atom_info(iatom,i)%coord(2)=mod(atom_info(iatom,i)%coord(2),boxsize(2))
        IF (atom_info(iatom,i)%coord(2) < 0) THEN
           atom_info(iatom,i)%coord(2) = atom_info(iatom,i)%coord(2) + boxsize(2) 
        ENDIF
        atom_info(iatom,i)%coord(3)=mod(atom_info(iatom,i)%coord(3),boxsize(3))
        IF (atom_info(iatom,i)%coord(3) < 0) THEN
           atom_info(iatom,i)%coord(3) = atom_info(iatom,i)%coord(3) + boxsize(3) 
        ENDIF

        WRITE(10,131) TRIM(atom_info(iatom, i)%atom_name), &
        atom_info(iatom,i)%coord(1), &
        atom_info(iatom,i)%coord(2), &
        atom_info(iatom,i)%coord(3)

    ENDDO

  ENDDO step

  !WRITE(*,*)"SAMPLE3 nb_divx = ", nb_divx

  write(6,*)'Sampled trajectory is written: ', sampled_pos_filename
  close(10)

END SUBROUTINE sample_and_recenter_format2
