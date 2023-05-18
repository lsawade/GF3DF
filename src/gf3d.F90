
module gf3d

  use constants, only: NCHANNELS, orientation, channels, MAX_STRING_LEN
  use ctypes, only: t_GF, t_source
  use calendar, only: julian_day
  use io, only: read_GF, print_GF, free_GF, write_output_SAC
  use utils, only: throwerror, scaleM
  use setup_source_location, only: setup_point_search_arrays
  use source_location, only: locate_sources
  use interpolation, only: interpolateMT
  use hdf5

  private
  public :: &
    read_GF, &
    print_GF, &
    t_GF, &
    free_GF, &
    throwerror, &
    setup_point_search_arrays, &
    t_source, &
    locate_sources, &
    scaleM, &
    julian_day, &
    interpolateMT, &
    write_output_SAC, &
    NCHANNELS, orientation, channels, MAX_STRING_LEN



end module gf3d