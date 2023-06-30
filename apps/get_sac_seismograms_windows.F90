program get_sac_seismograms_windows
!
!  Auto source file test
!

  use gf3d, only: write_seismograms_winsta, get_args, throwerror, init_log, finalize_log

  ! variable names
  character(len=256) :: hdf5_filename, source_filename, outputdir! input variable
  character(len=256), dimension(:), allocatable :: args
  character(len=256) :: network
  character(len=256) :: station
  double precision :: t0
  double precision :: dt
  integer :: it0
  integer :: itf
  logical :: OUTPUT_SEISMOS_SAC_BINARY = .true.
  logical :: OUTPUT_SEISMOS_SAC_ALPHANUM = .false.

  call init_log()

  args = get_args()

  ! Evaluate Args
  if (size(args) .ne. 9) then
    write (*,*) " "
    write (*,*) " "
    write (*,*) "Usage:"
    write (*,*) "----------------------------"
    write (*,*) " "
    write (*,*) "  gf3d-write-seismograms <H5file> <CMTSOLUTION> <output_dir> <network> <station> <t0> <dt> <it0> <itf>"
    write (*,*) " "
    write (*,*) "     H5file      = Subset file"
    write (*,*) "     CMTSOLUTION = CMTSOLUTION"
    write (*,*) "     output_dir  = OUTPUT_DIR"
    write (*,*) "     network     = network name"
    write (*,*) "     station     = station name "
    write (*,*) "     t0          = start time with respect to CMTSOLUTION"
    write (*,*) "     dt          = sampling interval"
    write (*,*) "     it0         = index of first sample"
    write (*,*) "     itf         = index of last sample"
    write (*,*) " "

    call throwerror(-1, "  *** gf3d-write-seismograms takes exactly 9 arguments.")

  endif


  ! which gives the filename
  hdf5_filename = args(1)
  source_filename = args(2)
  outputdir = args(3)

  network = args(4)
  station = args(5)
  read(args(6),*) t0
  read(args(7),*) dt
  read(args(8),*) it0
  read(args(9),*) itf

  write(*,*) "Args: ", trim(hdf5_filename), " ", trim(source_filename),  " ", &
                       trim(outputdir), " ", " ", trim(network), " ", trim(station), &
                       t0, dt, it0, itf

  ! Write seismograms
  call write_seismograms_winsta(&
      hdf5_filename, source_filename, outputdir, &
      network, station, &
      t0, dt, it0, itf, &
      OUTPUT_SEISMOS_SAC_ALPHANUM, &
      OUTPUT_SEISMOS_SAC_BINARY)

  call finalize_log()


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

  ! call cpu_time(finish)
  ! print '("Time = ",f6.3," seconds.")',finish-start
  ! print '("AvgTime = ",f6.3," seconds.")',(finish-start)/niter


 end program get_sac_seismograms_windows
