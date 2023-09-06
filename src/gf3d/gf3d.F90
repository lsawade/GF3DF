
module gf3d

  use gf, only: read_gf, print_GF, free_GF, t_GF
  use sources, only: read_cmt, print_source, t_source
  use utils, only: get_args, throwerror, nextpower2, init_log, finalize_log
  use stf, only: get_stf
  use sac, only: write_output_SAC
  use interpolation, only: spline1d, interp1d
  use seismograms, only: get_seismograms, write_seismograms
  use seismograms_windows, only: get_seismograms_winsta, write_seismograms_winsta

  ! use setup_source_location, only: setup_point_search_arrays
  use fftpack

  ! High level functions
  ! use gf3d_write_seismograms, only: write_seismograms
  ! use gf3d_get_seismograms, only: get_seismograms


  private
  public :: &
    t_GF, read_GF, print_GF, &
    t_source, read_cmt, print_source, &
    get_seismograms, &
    write_seismograms, &
    write_seismograms_winsta, &
    setup_point_search_arrays, &
    get_args, throwerror, nextpower2, &
    get_stf, &
    write_output_SAC, &
    interp1d, spline1d,  &
    init_log, finalize_log, &
    fftpack


end module gf3d