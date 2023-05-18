


program open_file

  use gf3d, only: read_GF, print_GF, t_GF, free_GF, throwerror, &
                  setup_point_search_arrays, t_source, locate_sources, &
                  scaleM, interpolateMT, write_output_SAC, &
                  NCHANNELS, orientation, channels, MAX_STRING_LEN, julian_day

  ! variable names
  character(len=65) :: filename ! input variable
  integer :: ix, num_args
  integer :: i, j, k, icomp
  integer(kind=8) :: iglob
  character(len=20), dimension(:), allocatable :: args
  type(t_source), dimension(1), target :: sources
  double precision, dimension(:,:,:), allocatable :: seismograms
  double precision, dimension(:,:,:,:,:,:,:), allocatable :: displacement
  character(len=MAX_STRING_LEN), parameter :: OUTPUT_DIR = './OUTPUT/'
  character(len=MAX_STRING_LEN) :: sisname
  character(len=MAX_STRING_LEN) :: model = "GLAD-M25"
  logical :: OUTPUT_SEISMOS_SAC_ALPHANUM = .true.
  logical :: OUTPUT_SEISMOS_SAC_BINARY = .true.

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
  sources(1)%Mrr = dble( 1.950000E+28)
  sources(1)%Mtt = dble(-4.360000E+26)
  sources(1)%Mpp = dble(-1.910000E+28)
  sources(1)%Mrt = dble( 7.420000E+27)
  sources(1)%Mrp = dble(-2.480000E+28)
  sources(1)%Mtp = dble( 9.420000E+26)
  sources(1)%eventname = "C201509162254A"
  sources(1)%time_shift = 49.9800
  sources(1)%hdur = 33.4000
  sources(1)%year = 2015
  sources(1)%month = 9
  sources(1)%day = 16
  sources(1)%jda = julian_day(sources(1)%year, sources(1)%month, sources(1)%day)
  sources(1)%hour = 22
  sources(1)%minute = 54
  sources(1)%second = 32.90

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

  allocate(seismograms(size(GF%displacement,1), 3, GF%nsteps))
  allocate(displacement(size(GF%displacement,1), 3, 3, GF%ngllx,GF%nglly,GF%ngllz,GF%nsteps))

  write(*,*) "SIZE  ", size(GF%displacement,1)
  ! Initialize arrays
  seismograms(:,:,:) = 0.d0
  displacement(:,:,:,:,:,:,:) = 0.d0

  !
  do k=1,GF%ngllz
    do j=1,GF%nglly
      do i=1,GF%ngllx
        iglob = GF%ibool(i,j,k, sources(1)%ispec)
        displacement(:,:,:,i,j,k, :) = GF%displacement(:,:,:,iglob,:)
      enddo
    enddo
  enddo

  call interpolateMT(&
    displacement, size(GF%displacement,1), &
    GF%nsteps, GF%ngllx, GF%nglly, GF%ngllz, GF%xigll, GF%yigll, GF%zigll, &
    sources(1)%Mxx, sources(1)%Myy, sources(1)%Mzz, sources(1)%Mxy, sources(1)%Mxz, sources(1)%Myz, &
    sources(1)%xi, sources(1)%eta,sources(1)%gamma, &
    sources(1)%xix, sources(1)%xiy, sources(1)%xiz, &
    sources(1)%etax, sources(1)%etay, sources(1)%etaz, &
    sources(1)%gammax, sources(1)%gammay, sources(1)%gammaz, &
    seismograms)



  do k=1,size(GF%displacement,1)
    write(*,*) trim(GF%networks(k)), ".", trim(GF%stations(k))

    do icomp=1,NCHANNELS

      write(sisname,"(a,'.',a)")
      write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
        trim(GF%networks(k)), trim(GF%stations(k)), channels(icomp)

      call write_output_SAC(&
        seismograms(k,icomp,:), &
        orientation(icomp), &
        sisname, &
        channels(icomp), &
        sources(1)%year, &
        sources(1)%jda, &
        sources(1)%hour, &
        sources(1)%minute, &
        sources(1)%second, &
        sources(1)%time_shift, &
        dble(0.0), &
        sources(1)%eventname, &
        sources(1)%latitude, &
        sources(1)%longitude, &
        sources(1)%depth, &
        sources(1)%hdur, &
        GF%stations(k), &
        GF%networks(k), &
        GF%latitudes(k), &
        GF%longitudes(k), &
        dble(0.0), &
        GF%burials(k), &
        GF%dt, &
        dble(0.0), &
        GF%nsteps, &
        OUTPUT_SEISMOS_SAC_ALPHANUM, &
        OUTPUT_SEISMOS_SAC_BINARY, &
        model, &
        OUTPUT_DIR)
    enddo
  enddo

  call free_GF(GF)
  deallocate(seismograms)

 end program open_file
