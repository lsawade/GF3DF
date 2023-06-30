program get_sac_seismograms
!
!  Auto source file test
!

  use gf3d, only: write_seismograms, get_args, throwerror, init_log, finalize_log

  ! variable names
  character(len=256) :: hdf5_filename, source_filename, outputdir! input variable
  character(len=256), dimension(:), allocatable :: args

  logical :: OUTPUT_SEISMOS_SAC_ALPHANUM = .false.
  logical :: OUTPUT_SEISMOS_SAC_BINARY = .false.

  call init_log()

  args = get_args()

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
      if (trim(args(4)) == "0") then
      OUTPUT_SEISMOS_SAC_BINARY = .true.
    else if (args(4) == "1") then
      OUTPUT_SEISMOS_SAC_ALPHANUM = .true.
    else if (args(4) == "2") then
      OUTPUT_SEISMOS_SAC_BINARY = .true.
      OUTPUT_SEISMOS_SAC_ALPHANUM = .true.
    else
      call throwerror(-1, "output type argument not recognized." )
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

  call finalize_log()

  ! real :: start, finish
  ! call cpu_time(start)

  !    put code to test here

  ! I haven't set the compilation flag for open mp yet...
  !$OMP PARALLEL DO
  ! do iter=1,niter


  !enddo
  !$OMP END PARALLEL DO

  ! call cpu_time(finish)
  ! print '("Time = ",f6.3," seconds.")',finish-start
  ! print '("AvgTime = ",f6.3," seconds.")',(finish-start)/niter


 end program get_sac_seismograms
