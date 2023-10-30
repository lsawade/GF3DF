


program open_file

  use gf3d, only: read_GF, print_GF, free_GF, throwerror, t_GF, get_args

  ! variable names
  character(len=65) :: filename ! input variable
  integer :: ix, num_args, k
  character(len=256), dimension(:), allocatable :: args

  type(t_GF) :: GF

  args = get_args()
  write (*,*) args
  if (size(args) /= 1) call throwerror(5, "*** gf3d takes exactly one argument which is an hdf5 file")

  ! which gives the filename
  filename = args(1)

  write(*,*) 'Reading Green Function file: ', trim(filename)

  ! Read Green Function file
  GF = read_GF(filename)

  ! Print header
  call print_GF(GF)

  ! call setup_point_search_arrays(GF)
  write (*,*) '# network, station, latitude, longitude, burial'
  do k=1,size(GF%displacement, 1)
    write(*,*) trim(GF%networks(k)),", ",trim(GF%stations(k)), ", ", GF%latitudes(k), ", ", GF%longitudes(k), ", ", GF%burials(k)
  enddo


  call free_GF(GF)


 end program open_file
