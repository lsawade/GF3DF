


program open_file

  use gf3d, only: read_GF, print_GF, t_GF, free_GF, throwerror

  ! variable names
  character(len=65) :: filename ! input variable
  integer :: ix, num_args
  character(len=20), dimension(:), allocatable :: args

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

  ! call setup_point_search_arrays(GF)

  call free_GF(GF)


 end program open_file
