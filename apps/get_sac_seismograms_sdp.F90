program get_sac_seismograms_sdp
!
!  Auto source file test
!

  use gf3d, only: write_seismograms, get_args, throwerror, init_log, finalize_log

  ! variable names
  character(len=256) :: hdf5_filename, source_filename, outputdir! input variable
  integer :: itypsokern
  character(len=256), dimension(:), allocatable :: args

  logical :: OUTPUT_SEISMOS_SAC_ALPHANUM = .false.
  logical :: OUTPUT_SEISMOS_SAC_BINARY = .false.

  call init_log()

  args = get_args()

  ! Evaluate Args
  if ((size(args) < 3) .or. (size(args) > 5)) then
    write (*,*) " "
    write (*,*) " Writes seismograms and partials to outputdirectory"
    write (*,*) " "
    write (*,*) "Usage:"
    write (*,*) "------"
    write (*,*) " "
    write (*,*) "  gf3d-write-seismograms <H5file> <CMTSOLUTION> <output_dir> <itypsokern> <outputtype>"
    write (*,*) " "
    write (*,*) "Description:"
    write (*,*) "------------"
    write (*,*) " "
    write (*,*) " itypsokern -- if 0, the subroutine only returns the synthetic "
    write (*,*) "                 seismogram in synt"
    write (*,*) "               if 1, the subroutine returns the synthetic seismogram "
    write (*,*) "                 in synt, and the 6 moment-tensor kernels in pdarray "
    write (*,*) "               if 2, the subroutine returns the synthetic seismogram "
    write (*,*) "                 in synt, the 6 moment-tensor kernels pdarray (index 1-6 )"
    write (*,*) "                 and 4 centroid kernels in pdarray (index 7-10 )The last "
    write (*,*) "                 argument is optional and takes values 0,1, or 2"
    write (*,*) " "
    write (*,*) " outputtype -- 0 -- outputs SAC BINARY"
    write (*,*) "               1 -- outputs SAC ALPHANUMERIC"
    write (*,*) "               2 -- both"
    write (*,*) " "
    write (*,*) " "
    call throwerror(-1, "  *** gf3d-write-seismograms takes 3 or 4 arguments")

  else if (size(args) == 3) then
    OUTPUT_SEISMOS_SAC_BINARY = .true.
    itypsokern = 0

  else if (size(args) == 4) then
    OUTPUT_SEISMOS_SAC_BINARY = .true.
    read(args(4), *) itypsokern

  else if (size(args) == 5) then
      if (trim(args(5)) == "0") then
      OUTPUT_SEISMOS_SAC_BINARY = .true.
    else if (args(5) == "1") then
      OUTPUT_SEISMOS_SAC_ALPHANUM = .true.
    else if (args(5) == "2") then
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
      hdf5_filename, source_filename, outputdir, itypsokern, &
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


 end program get_sac_seismograms_sdp
