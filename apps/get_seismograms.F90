program open_file
!
!  Auto source file test
!

  use gf3d, only:

                  read_GF, print_GF, t_GF, free_GF, throwerror, &
                  locate_sources, &
                  scaleM, interpolateMT, write_output_SAC, &
                  NCHANNELS, orientation, channels, MAX_STRING_LEN, julian_day, &
                  PI, read_cmt, print_source

  use stf, only: get_stf, stf_convolution, correct_hdur
  use fftpack, only: rk, fft, ifft, fftfreq
  use utils, only: nextpower2, get_args


  ! variable names
  character(len=65) :: hdf5_filename, source_filename ! input variable
  character(len=20), dimension(:), allocatable :: args
  integer :: i, j, k, icomp
  integer(kind=8) :: iglob
  type(t_source), dimension(1), target :: sources
  double precision, dimension(:,:,:), allocatable :: seismograms
  double precision, dimension(:,:,:,:,:,:,:), allocatable :: displacement
  character(len=MAX_STRING_LEN), parameter :: OUTPUT_DIR = './OUTPUT/'
  character(len=MAX_STRING_LEN) :: sisname
  character(len=MAX_STRING_LEN) :: model = "GLAD-M25"
  logical :: OUTPUT_SEISMOS_SAC_ALPHANUM = .true.
  logical :: OUTPUT_SEISMOS_SAC_BINARY = .true.

  double precision :: t0, tc, hdur_diff. shift
  integer :: NP2
  double precision, dimension(:), allocatable :: t, stf
  complex(kind=rk), dimension(:), allocatable :: cstf
  double precision, dimension(:,:,:), allocatable :: convolution


  type(t_GF) :: GF


  args = get_args()
  if (size(args) /= 2) call throwerror(5, "  *** write-seismograms takes exactly 2 arguments which is an hdf5 file, and a source file and output dir ***  ")


  ! which gives the filename
  hdf5_filename = args(1)
  source_filename = args(2)
  source_filename = args(2)

  !
  write_seismograms()

  ! Read Green Function file
  GF = read_GF(hdf5_filename)

  ! Print header
  call print_GF(GF)

  ! Read cmt solution
  sources = read_cmt(source_filename)

  write(*,*) "------------------------------------------------"
  call print_source(sources(1), 1)
  write(*,*) "------------------------------------------------"

  call locate_sources(GF, sources)

  call print_source(sources(1), 2)

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


  hdur_diff = correct_hdur(sources(1)%hdur, GF%hdur)

  t0 = 0.d0
  tc = 200.d0
  t(:) = t0 + ((/(I, I=1, GF%nsteps, 1)/)-1) * GF%dt

  call get_stf(t, tc, hdur_diff, int(GF%nsteps, kind=4), stf)

  allocate(convolution(size(GF%displacement,1), 3, GF%nsteps))
  convolution(:,:,:) = stf_convolution(seismograms, stf, GF%dt, tc)


  ! real :: start, finish
  ! call cpu_time(start)
  !       ! put code to test here


  ! I haven't set the compilation flag for open mp yet...
  !$OMP PARALLEL DO
  ! do iter=1,niter

  do k=1,size(GF%displacement,1)

    write(*,*) trim(GF%networks(k)), ".", trim(GF%stations(k))

    do icomp=1,NCHANNELS

      write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
      trim(GF%networks(k)), trim(GF%stations(k)), trim(channels(icomp))

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

  !enddo
  !$OMP END PARALLEL DO

  ! call cpu_time(finish)
  ! print '("Time = ",f6.3," seconds.")',finish-start
  ! print '("AvgTime = ",f6.3," seconds.")',(finish-start)/niter

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
  write(*,*) 'Nstat', size(GF%stations)
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
