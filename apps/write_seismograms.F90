program open_file

  use gf3d, only: read_GF, print_GF, t_GF, free_GF, throwerror, &
                  setup_point_search_arrays, t_source, locate_sources, &
                  scaleM, interpolateMT, write_output_SAC, &
                  NCHANNELS, orientation, channels, MAX_STRING_LEN, julian_day, &
                  PI, read_cmt, print_source

  use stf, only: get_stf, stf_convolution
  use fftpack, only: rk, fft, ifft, fftfreq
  use utils, only: nextpower2
  ! variable names
  character(len=65) :: hdf5_filename, source_filename ! input variable
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

  double precision :: t0, tc, hdur, hdur_conv, hdur_diff
  integer :: NP2
  double precision, dimension(:), allocatable :: t, stf, stf2
  complex(kind=rk), dimension(:), allocatable :: cstf
  double precision :: shift
  complex(kind=rk), dimension(:), allocatable :: pshift
  double precision, dimension(:), allocatable :: freqs
  complex(kind=rk), dimension(:), allocatable :: convo
  complex(kind=rk), dimension(:), allocatable :: cseis
  double precision, dimension(:), allocatable :: outseis
  double precision, dimension(:,:,:), allocatable :: convolution

  ! for the solution in time domain
  ! integer, parameter :: iratio = 32
  ! integer, parameter :: nfreq = 524288
  ! integer, parameter :: nt = iratio * nfreq

  ! integer :: it

  ! real :: wsave(4*nt+15)
  ! complex :: c(nt)



  type(t_GF) :: GF

  ! filename = '../single_element.h5' ! refer to line 5 of hdf5.js
  num_args = command_argument_count()
  if (num_args /= 2) call throwerror(5, "*** gf3d takes exactly 2 arguments which is an hdf5 file, and a source file")

  ! I've omitted checking the return status of the allocation
  allocate(args(num_args))

  ! Actually read each argument.
  do ix = 1, num_args
    call get_command_argument(ix,args(ix))
  end do

  ! which gives the filename
  hdf5_filename = args(1)
  source_filename = args(2)

  ! Read Green Function file
  GF = read_GF(hdf5_filename)

  ! Read cmt solution
  sources(:) = read_cmt(source_filename)

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
  sources(1)%force = .false.
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
  sources(1)%time_shift = 50.0000
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

  call print_source(sources(1), 1)

  allocate(t(GF%nsteps),stf(GF%nsteps), stf2(GF%nsteps))
  allocate(seismograms(size(GF%displacement,1), 3, GF%nsteps))
  allocate(displacement(size(GF%displacement,1), 3, 3, GF%ngllx,GF%nglly,GF%ngllz,GF%nsteps))

  write(*,*) "SIZE  ", size(GF%displacement,1)
  ! Initialize arrays
  seismograms(:,:,:) = 0.d0
  displacement(:,:,:,:,:,:,:) = 0.d0
  t(:) = 0.d0
  stf(:) = 0.d0

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



  if ((sources(1)%hdur / 1.628) ** 2 .le. GF%hdur**2) then

      hdur_diff = 0.000001

      write (*,*) &
          "Requested half duration smaller than what was simulated.\n", &
          "Half duration set to ", hdur_diff," s to simulate a Heaviside function."

      hdur_conv = sqrt(GF%hdur**2 - (sources(1)%hdur / 1.628)**2)

      write (*,*) "Try convolving your seismogram with a Gaussian with ", &
                  hdur_conv, " standard deviation."
  else
    hdur_diff = sqrt((sources(1)%hdur / 1.628)**2 - GF%hdur**2)
  endif


  t0 = 0.d0
  tc = 200.d0
  t(:) = t0 + ((/(I, I=1, GF%nsteps, 1)/)-1) * GF%dt
  ! do i=1,GF%nsteps
  !   t(i) = t0 + (i-1)*GF%dt
  ! enddo

  call get_stf(t, tc, hdur_diff, int(GF%nsteps, kind=4), stf)


  shift = -200.0
  allocate(convolution(size(GF%displacement,1), 3, GF%nsteps))
  convolution(:,:,:) = stf_convolution(seismograms, stf2, GF%dt, shift)


  ! real :: start, finish
  ! call cpu_time(start)
  !       ! put code to test here

  ! do iter=1,niter

  do k=1,size(GF%displacement,1)
    ! write(*,*) trim(GF%networks(k)), ".", trim(GF%stations(k))


    write(sisname,"(a,'.',a)")
    write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
      trim(GF%networks(k)), trim(GF%stations(k)), channels(icomp)

    do icomp=1,NCHANNELS

      call write_output_SAC(&
        ! seismograms(k,icomp,:), &
        convolution(k,icomp,:), &
        orientation(icomp), &
        sisname, &
        channels(icomp), &
        sources(1)%year, &
        sources(1)%jda, &
        sources(1)%hour, &
        sources(1)%minute, &
        sources(1)%second, &
        sources(1)%time_shift, &
        dble(GF%tc), &
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


  call exit()


  call get_stf(t, tc, hdur_diff, int(GF%nsteps, kind=4), stf)

  write(sisname,"(a,'.',a)")
  write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
    "STF", "ERF", "TIM"


  call write_output_SAC(&
        stf, &
        orientation(1), &
        sisname, &
        channels(1), &
        sources(1)%year, &
        sources(1)%jda, &
        sources(1)%hour, &
        sources(1)%minute, &
        sources(1)%second, &
        sources(1)%time_shift, &
        dble(GF%tc), &
        sources(1)%eventname, &
        sources(1)%latitude, &
        sources(1)%longitude, &
        sources(1)%depth, &
        sources(1)%hdur, &
        GF%stations(1), &
        GF%networks(1), &
        GF%latitudes(1), &
        GF%longitudes(1), &
        dble(0.0), &
        GF%burials(1), &
        GF%dt, &
        dble(0.0), &
        GF%nsteps, &
        OUTPUT_SEISMOS_SAC_ALPHANUM, &
        OUTPUT_SEISMOS_SAC_BINARY, &
        model, &
        OUTPUT_DIR)


  NP2 = nextpower2(2 * int(GF%nsteps, kind=4))

  cstf = fft(cmplx(stf,kind=rk), NP2)

  write(*,*) 'Sizes', GF%nsteps, NP2

  stf = dble(abs(cstf))

  write(sisname,"(a,'.',a)")
  write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
    "STF", "ERF", "FRE"

  call write_output_SAC(&
        stf, &
        orientation(1), &
        sisname, &
        channels(1), &
        sources(1)%year, &
        sources(1)%jda, &
        sources(1)%hour, &
        sources(1)%minute, &
        sources(1)%second, &
        sources(1)%time_shift, &
        dble(GF%tc), &
        sources(1)%eventname, &
        sources(1)%latitude, &
        sources(1)%longitude, &
        sources(1)%depth, &
        sources(1)%hdur, &
        GF%stations(1), &
        GF%networks(1), &
        GF%latitudes(1), &
        GF%longitudes(1), &
        dble(0.0), &
        GF%burials(1), &
        GF%dt, &
        dble(0.0), &
        GF%nsteps, &
        OUTPUT_SEISMOS_SAC_ALPHANUM, &
        OUTPUT_SEISMOS_SAC_BINARY, &
        model, &
        OUTPUT_DIR)


  cstf = ifft(cstf, NP2)

  stf = dble(real(cstf, kind=8))

  write(sisname,"(a,'.',a)")
  write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
    "STF", "ERF", "IFR"

  call write_output_SAC(&
        stf, &
        orientation(1), &
        sisname, &
        channels(1), &
        sources(1)%year, &
        sources(1)%jda, &
        sources(1)%hour, &
        sources(1)%minute, &
        sources(1)%second, &
        sources(1)%time_shift, &
        dble(GF%tc), &
        sources(1)%eventname, &
        sources(1)%latitude, &
        sources(1)%longitude, &
        sources(1)%depth, &
        sources(1)%hdur, &
        GF%stations(1), &
        GF%networks(1), &
        GF%latitudes(1), &
        GF%longitudes(1), &
        dble(0.0), &
        GF%burials(1), &
        GF%dt, &
        dble(0.0), &
        GF%nsteps, &
        OUTPUT_SEISMOS_SAC_ALPHANUM, &
        OUTPUT_SEISMOS_SAC_BINARY, &
        model, &
        OUTPUT_DIR)



  ! call cpu_time(finish)
  ! print '("Time = ",f6.3," seconds.")',finish-start
  ! print '("AvgTime = ",f6.3," seconds.")',(finish-start)/niter

  call free_GF(GF)
  deallocate(seismograms)

 end program open_file
