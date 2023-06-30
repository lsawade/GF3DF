
module gf3d

  use gf, only: read_gf, print_GF, free_GF
  use sources, only: read_cmt, print_source, t_source
  use utils, only: get_args, throwerror, nextpower2, init_log, finalize_log
  use stf, only: get_stf
  use sac, only: write_output_SAC
  use interpolation, only: spline1d, interp1d
  use fftpack

  ! High level functions
  ! use gf3d_write_seismograms, only: write_seismograms
  ! use gf3d_get_seismograms, only: get_seismograms


  private
  public :: &
    read_GF, read_cmt, &
    print_GF, print_source, t_source, &
    write_seismograms, &
    write_seismograms_winsta, &
    get_args, throwerror, nextpower2, &
    get_stf, &
    write_output_SAC, &
    interp1d, spline1d,  &
    init_log, finalize_log, &
    fftpack

  interface get_seismograms
    module subroutine get_seismograms(GF, sources, superseismograms)
      use gf, only: t_GF
      use sources, only: t_source
      type(t_GF), intent(in) :: GF
      type(t_source), dimension(:), intent(inout) :: sources
      double precision, dimension(:,:,:), allocatable, intent(out) :: superseismograms
    end subroutine
  end interface get_seismograms

  interface write_seismograms
    module subroutine write_seismograms(&
      GF_filename, source_filename, output_dir, &
      OUTPUT_SEISMOS_SAC_ALPHANUM, &
      OUTPUT_SEISMOS_SAC_BINARY)
      character(len=*), intent(in) :: GF_filename, source_filename, output_dir
      logical, intent(in) :: OUTPUT_SEISMOS_SAC_ALPHANUM
      logical, intent(in) :: OUTPUT_SEISMOS_SAC_BINARY
    end subroutine
  end interface write_seismograms

  interface get_seismograms_winsta
    module subroutine get_seismograms_winsta(&
      GF, sources, &
      istat, &
      t0, dt, it0, itf,  &
      superseismograms)
      use gf, only: t_GF
      use sources, only: t_source

      type(t_GF), intent(in) :: GF
      type(t_source), dimension(:), intent(inout) :: sources
      integer, intent(in) :: istat
      double precision, intent(in) :: t0
      double precision, intent(in) :: dt
      integer, intent(in) :: it0
      integer, intent(in) :: itf
      double precision, dimension(:,:,:), allocatable, intent(out) :: superseismograms
    end subroutine
  end interface get_seismograms_winsta

  interface write_seismograms_winsta
    module subroutine write_seismograms_winsta(&
      GF_filename, source_filename, output_dir, &
      network, station, &
      t0, dt, it0, itf,  &
      OUTPUT_SEISMOS_SAC_ALPHANUM, &
      OUTPUT_SEISMOS_SAC_BINARY)
      character(len=*), intent(in) :: GF_filename, source_filename, output_dir
      character(len=*), intent(in) :: network, station
      double precision, intent(in) :: t0
      double precision, intent(in) :: dt
      integer, intent(in) :: it0
      integer, intent(in) :: itf
      logical, intent(in) :: OUTPUT_SEISMOS_SAC_ALPHANUM
      logical, intent(in) :: OUTPUT_SEISMOS_SAC_BINARY
    end subroutine
  end interface write_seismograms_winsta




end module gf3d