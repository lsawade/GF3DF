program print_cmt

  use sources, only: t_source, read_cmt, print_source
  use utils, only: get_args, throwerror

  ! Types
  integer :: i
  type(t_source), dimension(:), allocatable :: sources
  character(len=20), dimension(:), allocatable :: args

  args = get_args()

  if (size(args) .ne. 1) call throwerror(5, "Only takes one argument")

  ! Read file
  sources = read_cmt(trim(args(1)))

  ! Print source
  do i=1,size(sources)
    call print_source(sources(i), 1)
  enddo


end program print_cmt