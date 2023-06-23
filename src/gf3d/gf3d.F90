
module gf3d

  use gf, only: read_gf, print_GF, free_GF
  use sources, only: read_cmt, print_source, t_source
  use utils, only: get_args, throwerror
  use gf3d_write_seismograms, only: write_seismograms
  use gf3d_get_seismograms, only: get_seismograms

  private
  public :: &
    read_GF, read_cmt, &
    print_GF, print_source, t_source, &
    write_seismograms, &
    get_args, throwerror

end module gf3d