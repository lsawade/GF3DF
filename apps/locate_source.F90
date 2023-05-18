


program open_file

  use gf3d, only: read_GF, print_GF, t_GF, free_GF, throwerror, &
                  setup_point_search_arrays, t_source, locate_sources, &
                  scaleM

  ! variable names
  character(len=65) :: filename ! input variable
  integer :: ix, num_args
  character(len=20), dimension(:), allocatable :: args
  type(t_source), dimension(1), target :: sources

  type(t_GF) :: GF

  ! filename = '../single_element.h5' ! refer to line 5 of hdf5.js
  num_args = command_argument_count()
  if (num_args /= 1) call throwerror(5, "*** gf3d takes exactly one argument which is an hdf5 file")

  ! I've omitted checking the return status of the allocation
  allocate(args(num_args))

  ! Actually read each argument.
  do ix = 1, num_args
    call get_command_argument(ix,args(ix))
  end do

  ! which gives the filename
  filename = args(1)

  ! Read Green Function file
  GF = read_GF(filename)

  ! Print header
  call print_GF(GF)

  call setup_point_search_arrays(GF)



  ! ---------------------------------------------------------------------------
  ! Using this source
  ! PDEW2015  9 16 22 54 32.90 -31.5700  -71.6700  22.4 0.0 8.3 NEAR COAST OF CENTRAL CH
  ! event name:     201509162254A
  ! time shift:     49.9800
  ! half duration:  33.4000
  ! latitude:      -31.1300
  ! longitude:     -72.0900
  ! depth:          17.3500
  ! Mrr:       1.950000e+28
  ! Mtt:      -4.360000e+26
  ! Mpp:      -1.910000e+28
  ! Mrt:       7.420000e+27
  ! Mrp:      -2.480000e+28
  ! Mtp:       9.420000e+26
  ! ---------------------------------------------------------------------------

  ! Read all the sources
  ! call read_source_locations(sources)

  ! Define Source
  sources(1)%latitude = -31.1300
  sources(1)%longitude = -72.0900
  sources(1)%depth = 17.3500
  sources(1)%Mrr = scaleM(dble( 1.950000E+28))
  sources(1)%Mtt = scaleM(dble(-4.360000E+26))
  sources(1)%Mpp = scaleM(dble(-1.910000E+28))
  sources(1)%Mrt = scaleM(dble( 7.420000E+27))
  sources(1)%Mrp = scaleM(dble(-2.480000E+28))
  sources(1)%Mtp = scaleM(dble( 9.420000E+26))
  sources(1)%time_shift = 49.9800
  sources(1)%hdur = 33.4000

  write(*,*) "GF%bool"
  write(*,*) GF%ibool

  call locate_sources(GF, sources)

  write(*,*) "Located source values"
  write(*,*) "------------------------------------------------"
  write(*,*) "Mxx    ", sources(1)%Mxx
  write(*,*) "Myy    ", sources(1)%Myy
  write(*,*) "Mzz    ", sources(1)%Mzz
  write(*,*) "Mxy    ", sources(1)%Mxy
  write(*,*) "Mxz    ", sources(1)%Mxz
  write(*,*) "Myz    ", sources(1)%Myz
  write(*,*) "x      ", sources(1)%x
  write(*,*) "y      ", sources(1)%y
  write(*,*) "z      ", sources(1)%z
  write(*,*) "xix    ", sources(1)%xix
  write(*,*) "xiy    ", sources(1)%xiy
  write(*,*) "xiz    ", sources(1)%xiz
  write(*,*) "etax   ", sources(1)%etax
  write(*,*) "etay   ", sources(1)%etay
  write(*,*) "etaz   ", sources(1)%etaz
  write(*,*) "gammax ", sources(1)%gammax
  write(*,*) "gammay ", sources(1)%gammay
  write(*,*) "gammaz ", sources(1)%gammaz
  write(*,*) "------------------------------------------------"

  call free_GF(GF)


 end program open_file
