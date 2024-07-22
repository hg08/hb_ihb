!The aim of this program is to have information about density in a cell
!Three information are provided: -The density in g.m-3
!                                -The density in at.A-1
!                                -The integral of at.A-1
program main
  implicit none

  !Definition of a structure to store the informations about common elements (already used atoms) and new elements (to define when the program runs)                  
  type info_atoms
     character(len=3) name
     double precision mw
  end type info_atoms


  integer step,stepi,stepf,naxis,nstep,ndiv,at,nat,nattot,io,i,j,kindmax
  integer, dimension(:), allocatable :: tabsel
  double precision, parameter :: Na = 6.02214076*10.0**(-7) !-7 is writen instead of 23 because it will devided by m**-3.Ang**-3 (10**-10)**3 
  double precision, parameter :: rate_m2cm = 10**6   ! 1 m^3 = 10^6 cm^3
  double precision volume
  double precision, dimension(3):: pos,paramcell
  double precision, dimension(:), allocatable :: tabdens,tabat,tabsum,tabweight
  character char
  character(len=3),dimension(:),allocatable::tabsymbol
  character(len=200)::filein,fileout,symbol,selection,chartmp
  type (info_atoms), dimension(:), allocatable :: tab_info

  filein=""
  fileout=""

  !==============================
  !Information about common atoms
  !==============================
  kindmax=10
  allocate (tab_info(kindmax))
  tab_info(1)%name='H'  ; tab_info(1)%mw=1.0d0  
  tab_info(2)%name='Li' ; tab_info(2)%mw=7.0d0  
  tab_info(3)%name='C'  ; tab_info(3)%mw=12.0d0 
  tab_info(4)%name='N'  ; tab_info(4)%mw=14.0d0 
  tab_info(5)%name='O'  ; tab_info(5)%mw=16.0d0 
  tab_info(6)%name='F'  ; tab_info(6)%mw=19.0d0 
  tab_info(7)%name='Cl' ; tab_info(7)%mw=35.453d0
  tab_info(8)%name='Ca' ; tab_info(8)%mw=40.1d0 
  tab_info(9)%name='Na' ; tab_info(9)%mw=22.99d0
  tab_info(10)%name='I' ; tab_info(10)%mw=126.9d0


  !=================
  !Basic information
  !=================

  write(*,*)"Name of the input file:"!pos-file
  read(*,*)filein
  open(10,file=filein,iostat=io)
  if(io<0)then
     write(*,*) "File does not exist."
     stop
  endif

  read(10,*)nattot
  write(*,*)"There are ",nattot," atoms."
  allocate(tabsel(nattot),tabweight(nattot))
  read(10,*)char,char, step
  write(*,*)"The starting step is ",step,"."

  write(*,*)"Name of the output file:"
  read(*,'(a)')fileout
  open(11,file=fileout)
  write(*,*)"Number of division for the density (default 500):"
  read(*,'(i30)')ndiv
  if(ndiv==0)ndiv=501
  write(*,*)"You want ",ndiv," divisions."
  allocate(tabdens(ndiv),tabat(ndiv),tabsum(ndiv))
  tabdens=0
  tabat=0

  write(*,*)"Initial step (default",step,"):"
  read(*,'(i30)')stepi
  if((stepi==0).or.(stepi<step))stepi=step
  write(*,*)"You choose", stepi,"."

  write(*,*)"Final step (default 'end of file'):"
  read(*,'(i30)')stepf
  write(*,*)"You choose", stepf,"."
  write(*,*)

  write(*,*)"Cell parameters (Angstrom):"
  read(*,*)paramcell(1),paramcell(2),paramcell(3)
  write(*,*)"The length is ",paramcell(1),",",paramcell(2),",",paramcell(3),"."
  write(*,*)

  write(*,*)"Direction to study (x/y/z, default z):"
  read(*,'(a)')char
  if(char=="")char="z"
  write(*,*)"You chose the direction ",char,"."
  if(char=="x")naxis=1
  if(char=="y")naxis=2
  if(char=="z")naxis=3
  volume=paramcell(1)*paramcell(2)*paramcell(3)
  
  write(*,*)"Which atoms have to be selected? (ALL/NUMBER/NAME, default: ALL)"
  read(*,'(a)')selection
  if(selection=="")selection="ALL"
  write(*,*)"You have selected ",trim(selection),"."
  tabsel=0


  if(selection=="ALL")then
     tabsel=1!tabsel is a table which records the atoms which will be studied: 0=no studied, 1=studied
  elseif(selection=="NUMBER")then
     write(*,*)"How many atoms will be selected?"
     read(*,*)nat
     write(*,*)"Write the atoms number separated by a coma:"
     do i=1,nat-1
        read(*,'(i10)',advance='no')at
        tabsel(at)=1
     enddo
     read(*,'(i10)')at
     tabsel(at)=1
  elseif(selection=="NAME")then
     write(*,*)"How many kind of atoms will be selected?"
     read(*,*)nat
     allocate(tabsymbol(nat))
     write(*,*)"Write the symbol of the atoms:"
     !For the NAME keyword,tabsel is filled during the reading of the first step 
     !Only way that I found to separate each symbol (1, 2 or 3 characters) and to fill tabsymbol.
     read(*,'(a)')chartmp   
     do i=1,nat
        chartmp=adjustl(chartmp)
        tabsymbol(i)=chartmp(:index(chartmp," "))
        chartmp=chartmp(index(chartmp," "):)
     enddo
  endif

  !======================
  !Reading the first step
  !======================
  !1)tabsel is not filled in the case of selection=NAME but it is already done for ALL and NUMBER
  !2)tabweight is filled with the molecular weight of each atom
  !So we get these values thanks to the first step

  do i=1,nattot!Loop for the atoms
     read(10,*)symbol

     !1) if selection="NAME", we fill tabsel
     if(selection=="NAME")then
        do j=1,nat
           if(symbol==tabsymbol(j))tabsel(i)=1
        enddo
     endif

     !2)Filling tabweight
     do j=1,kindmax
        if(tab_info(j)%name==trim(symbol))then
           tabweight(i)=tab_info(j)%mw
           exit
        endif
     enddo
     if(j==kindmax+1)then
        write(*,*)"Be careful the weight of",symbol," is unknown!"
        exit
     endif

  enddo


  !Come back to the beginning of the file if necessary
  !And skipping the 2 header lines
  if(step==stepi)then
     REWIND(10)
  endif
  read(10,*)
  read(10,*)char,char, step

  !==================
  !Input file reading
  !==================
  nstep=0
  !Passing the uselles steps
  do while(step<stepi)
     do i=1,nattot
        read(10,*)
     enddo
     read(10,*)
     read(10,*)char,char,step
  enddo

  do while ((step<=stepf).or.(stepf==0))!Loop for the steps
     nstep=nstep+1

     do i=1,nattot!Loop for the atoms
        read(10,*)symbol,pos(1),pos(2),pos(3)

        !==============
        !Filing tabdens
        !==============        
        if(tabsel(i)==1)then
           j=nint( (pos(naxis)/paramcell(naxis) - floor(pos(naxis)/paramcell(naxis)) ) * (ndiv-1) ) + 1 !j=[1;ndiv], for density , each step we need to recenter (because of the table allocation).
           tabdens(j)=tabdens(j)+tabweight(i)
           tabat(j)=tabat(j)+1
        endif
     enddo

     read(10,*,iostat=io)
     if(io<0) exit!end if we reach the end of the input file
     read(10,*)char,char,step
     if(mod(step,5000)==0)write(*,*)"Step",step," recorded."

  enddo

  !============
  !Writing file
  !============
  write(*,*) "ndiv-1=",ndiv-1
  write(*,*) "nstep=",nstep
  write(*,*) "volume=",volume

  !Header
  write(11,*)"#Length(Ang)    Density(g.cm-3)      Density(nbat.Ang-1)    Integral(g.m-3)"

  !Body
  tabdens=tabdens*(ndiv-1)/nstep/Na/volume!Density in g.m-3
  !tabdens=tabdens /nstep/Na / (volume/(ndiv-1))!Density in g.m-3
  tabat=tabat*(ndiv-1)/nstep/paramcell(naxis)!Density in at.A-1
  tabsum(1)=tabat(1)!Integral of the density in at.A-1
  do i=2,ndiv
     tabsum(i)=tabsum(i-1)+tabat(i)
  enddo
  
  ! Convert to g/cm^3
  tabdens = tabdens/rate_m2cm

  do i=1,ndiv
     !Density
     write(11,'(f20.10,f20.10,f20.10,f20.10)')(i-1.)/(ndiv-1.)*paramcell(naxis),tabdens(i),tabat(i),tabsum(i)
  enddo


  deallocate(tabsel,tab_info,tabdens,tabweight,tabat,tabsum)
  if(selection=="NAME")deallocate(tabsymbol)

  close(10)
  close(11)


end program main
