program write_sac_seismograms
!
!  Auto source file test
!

  use gf3d, only: write_seismograms, get_args, throwerror

  ! variable names
  character(len=65) :: hdf5_filename, source_filename, outputdir! input variable
  character(len=20), dimension(:), allocatable :: args

  logical :: OUTPUT_SEISMOS_SAC_ALPHANUM = .false.
  logical :: OUTPUT_SEISMOS_SAC_BINARY = .false.

  args = get_args()

  write (*,*) args
  write (*,*) "size", size(args)

  ! Evaluate Args
  if ((size(args) < 3) .or. (size(args) > 4)) then
    write (*,*) " "
    write (*,*) " "
    write (*,*) "Usage:"
    write (*,*) "----------------------------"
    write (*,*) " "
    write (*,*) "  gf3d-write-seismograms <H5file> <CMTSOLUTION> <output_dir> [{0,1,2}]"
    write (*,*) " "
    write (*,*) "  The last argument is optional and takes values 0,1, or 2"
    write (*,*) " "
    write (*,*) "  0 -- outputs SAC BINARY"
    write (*,*) "  1 -- outputs SAC ALPHANUMERIC"
    write (*,*) "  2 -- both"
    write (*,*) " "
    write (*,*) " "
    call throwerror(-1, "  *** gf3d-write-seismograms takes 3 or 4 arguments")

  else if (size(args) == 3) then
    OUTPUT_SEISMOS_SAC_BINARY = .true.

  else if (size(args) == 4) then

    if (args(3) == "0") then
      OUTPUT_SEISMOS_SAC_BINARY = .true.
    else if (args(3) == "1") then
      OUTPUT_SEISMOS_SAC_ALPHANUM = .true.
    else if (args(3) == "2") then
      OUTPUT_SEISMOS_SAC_BINARY = .true.
      OUTPUT_SEISMOS_SAC_ALPHANUM = .true.
    endif

  endif


  ! which gives the filename
  hdf5_filename = args(1)
  source_filename = args(2)
  outputdir = args(3)

  ! Write seismograms
  call write_seismograms(&
      hdf5_filename, source_filename, outputdir, &
      OUTPUT_SEISMOS_SAC_ALPHANUM, OUTPUT_SEISMOS_SAC_BINARY)


  ! real :: start, finish
  ! call cpu_time(start)
  !       ! put code to test here


  ! I haven't set the compilation flag for open mp yet...
  !$OMP PARALLEL DO
  ! do iter=1,niter


  !enddo
  !$OMP END PARALLEL DO

  ! call cpu_time(finish)
  ! print '("Time = ",f6.3," seconds.")',finish-start
  ! print '("AvgTime = ",f6.3," seconds.")',(finish-start)/niter

  ! call get_stf(t, tc, hdur_diff, int(GF%nsteps, kind=4), stf)

  ! write(sisname,"(a,'.',a)")
  ! write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
  !   "STF", "ERF", "TIM"


  ! call write_output_SAC(&
  !       stf, &
  !       orientation(1), &
  !       sisname, &
  !       channels(1), &
  !       sources(1)%year, &
  !       sources(1)%jda, &
  !       sources(1)%hour, &
  !       sources(1)%minute, &
  !       sources(1)%second, &
  !       sources(1)%time_shift, &
  !       dble(GF%tc), &
  !       sources(1)%eventname, &
  !       sources(1)%latitude, &
  !       sources(1)%longitude, &
  !       sources(1)%depth, &
  !       sources(1)%hdur, &
  !       GF%stations(1), &
  !       GF%networks(1), &
  !       GF%latitudes(1), &
  !       GF%longitudes(1), &
  !       dble(0.0), &
  !       GF%burials(1), &
  !       GF%dt, &
  !       dble(0.0), &
  !       GF%nsteps, &
  !       OUTPUT_SEISMOS_SAC_ALPHANUM, &
  !       OUTPUT_SEISMOS_SAC_BINARY, &
  !       model, &
  !       OUTPUT_DIR)


  ! NP2 = nextpower2(2 * int(GF%nsteps, kind=4))

  ! cstf = fft(cmplx(stf,kind=rk), NP2)

  ! write(*,*) 'Sizes', GF%nsteps, NP2
  ! write(*,*) 'Nstat', size(GF%stations)
  ! stf = dble(abs(cstf))

  ! write(sisname,"(a,'.',a)")
  ! write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
  !   "STF", "ERF", "FRE"

  ! call write_output_SAC(&
  !       stf, &
  !       orientation(1), &
  !       sisname, &
  !       channels(1), &
  !       sources(1)%year, &
  !       sources(1)%jda, &
  !       sources(1)%hour, &
  !       sources(1)%minute, &
  !       sources(1)%second, &
  !       sources(1)%time_shift, &
  !       dble(GF%tc), &
  !       sources(1)%eventname, &
  !       sources(1)%latitude, &
  !       sources(1)%longitude, &
  !       sources(1)%depth, &
  !       sources(1)%hdur, &
  !       GF%stations(1), &
  !       GF%networks(1), &
  !       GF%latitudes(1), &
  !       GF%longitudes(1), &
  !       dble(0.0), &
  !       GF%burials(1), &
  !       GF%dt, &
  !       dble(0.0), &
  !       GF%nsteps, &
  !       OUTPUT_SEISMOS_SAC_ALPHANUM, &
  !       OUTPUT_SEISMOS_SAC_BINARY, &
  !       model, &
  !       OUTPUT_DIR)


  ! cstf = ifft(cstf, NP2)

  ! stf = dble(real(cstf, kind=8))

  ! write(sisname,"(a,'.',a)")
  ! write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
  !   "STF", "ERF", "IFR"

  ! call write_output_SAC(&
  !       stf, &
  !       orientation(1), &
  !       sisname, &
  !       channels(1), &
  !       sources(1)%year, &
  !       sources(1)%jda, &
  !       sources(1)%hour, &
  !       sources(1)%minute, &
  !       sources(1)%second, &
  !       sources(1)%time_shift, &
  !       dble(GF%tc), &
  !       sources(1)%eventname, &
  !       sources(1)%latitude, &
  !       sources(1)%longitude, &
  !       sources(1)%depth, &
  !       sources(1)%hdur, &
  !       GF%stations(1), &
  !       GF%networks(1), &
  !       GF%latitudes(1), &
  !       GF%longitudes(1), &
  !       dble(0.0), &
  !       GF%burials(1), &
  !       GF%dt, &
  !       dble(0.0), &
  !       GF%nsteps, &
  !       OUTPUT_SEISMOS_SAC_ALPHANUM, &
  !       OUTPUT_SEISMOS_SAC_BINARY, &
  !       model, &
  !       OUTPUT_DIR)



  ! call cpu_time(finish)
  ! print '("Time = ",f6.3," seconds.")',finish-start
  ! print '("AvgTime = ",f6.3," seconds.")',(finish-start)/niter


 end program write_sac_seismograms
