program get_sac_seismograms_sdp
!
!  Auto source file test
!

  use gf3d, only: get_seismograms, get_args, throwerror, init_log, finalize_log, &
                  t_source, t_GF, read_GF, setup_point_search_arrays

  ! Command line args
  character(len=256) :: hdf5_filename
  integer :: itypsokern, niter = 10
  character(len=256), dimension(:), allocatable :: args

  ! Important parameters
  type(t_GF):: GF
  type(t_source), dimension(1) :: sources
  double precision, dimension(:,:,:), allocatable :: synt
  double precision, dimension(:,:,:,:), allocatable :: dp

  ! Local things
  real :: start, finish

  call init_log()

  args = get_args()

  ! Evaluate Args
  if ((size(args) < 1) .or. (size(args) > 2)) then
    write (*,*) " "
    write (*,*) " Writes seismograms and partials to outputdirectory"
    write (*,*) " "
    write (*,*) "Usage:"
    write (*,*) "------"
    write (*,*) " "
    write (*,*) "  gf3d-write-seismograms <H5file> [<itypsokern>]"
    write (*,*) " "
    write (*,*) "Description:"
    write (*,*) "------------"
    write (*,*) " "
    write (*,*) " itypsokern -- if 0, the subroutine only returns the synthetic "
    write (*,*) "                 seismogram in synt. Default if omitted"
    write (*,*) "               if 1, the subroutine returns the synthetic seismogram "
    write (*,*) "                 in synt, and the 6 moment-tensor kernels in pdarray "
    write (*,*) "               if 2, the subroutine returns the synthetic seismogram "
    write (*,*) "                 in synt, the 6 moment-tensor kernels pdarray (index 1-6 )"
    write (*,*) "                 and 4 centroid kernels in pdarray (index 7-10 )The last "
    write (*,*) "                 argument is optional and takes values 0,1, or 2"
    write (*,*) " "
    write (*,*) " "
    call throwerror(-1, "  *** gf3d-write-seismograms takes 3 or 4 arguments")
  else if (size(args) == 1) then
    itypsokern = 0
  else if (size(args) == 2) then
    read(args(2), *) itypsokern
  endif

  ! Which gives the filename
  hdf5_filename = args(1)

  ! Given the filename we can first load the Green function
  GF = read_GF(hdf5_filename)

  ! Setup KDTree, this has to be done like this du to module structure
  ! I might change this later, but only if I have a better idea.
  call setup_point_search_arrays(GF)

  ! Define Source parameters (I assume here you have your own parameters)
  sources(1)%latitude = -31.1300
  sources(1)%longitude = -72.0900
  sources(1)%depth = 17.3500
  sources(1)%Mrr = dble( 1.950000E+28)
  sources(1)%Mtt = dble(-4.360000E+26)
  sources(1)%Mpp = dble(-1.910000E+28)
  sources(1)%Mrt = dble( 7.420000E+27)
  sources(1)%Mrp = dble(-2.480000E+28)
  sources(1)%Mtp = dble( 9.420000E+26)
  sources(1)%time_shift = 49.9800
  sources(1)%hdur = 33.4000

  ! Get start timestamp to compute average extraction time
  write(*,*) "Starting extraction of ", niter, " kernel sets."
  call cpu_time(start)

  ! Extract seismograms niter times
  do i=1, niter
    call get_seismograms(GF, sources, synt, dp, itypsokern)
    write (*,*) "  ... done ", i,"."
  enddo

  ! Get finish timestamp to compute average extraction time
  call cpu_time(finish)

  ! Print finalizing statement
  write(*,*) "Average time for extraction of a single source with itypsokern=", &
             itypsokern, "is  ", (finish-start)/niter,"s."

  ! Finish logging
  call finalize_log()

 end program get_sac_seismograms_sdp
