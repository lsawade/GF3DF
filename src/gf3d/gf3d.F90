
module gf3d

  use gf, only: read_gf, print_GF, free_GF, t_GF
  use sources, only: read_cmt, print_source, t_source
  use utils, only: get_args, throwerror, nextpower2, init_log, finalize_log
  use stf, only: get_stf
  use sac, only: write_output_SAC
  use interpolation, only: spline1d, interp1d
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
    ! setup_point_search_arrays, &
    get_args, throwerror, nextpower2, &
    get_stf, &
    write_output_SAC, &
    interp1d, spline1d,  &
    init_log, finalize_log
    ! fftpack

  ! interface interpolate_source
  !   module subroutine interpolate_source(GF, source, seismograms)
  !     use gf, only: t_GF
  !     use sources, only: t_source
  !     use interpolation, only: interpolateMT
  !     type(t_GF), intent(in) :: GF
  !     type(t_source), intent(in) :: source
  !     double precision, dimension(:,:,:) :: seismograms
  !   end subroutine interpolate_source
  ! end interface interpolate_source

  ! Get seismograms interface
  interface get_seismograms

    ! Simple get_seismograms usage
    module subroutine get_seismograms(GF, sources, superseismograms)
      use gf, only: t_GF
      use sources, only: t_source
      type(t_GF), intent(in) :: GF
      type(t_source), dimension(:), intent(inout) :: sources
      double precision, dimension(:,:,:), allocatable, intent(out) :: superseismograms
    end subroutine

    ! Get sdp as used in GCMT

    module subroutine get_sdp(GF, sources, synt, dp, itypsokern)
      use gf, only: t_GF
      use sources, only: t_source
      use stf, only: get_stf, stf_convolution, correct_hdur
      use interpolation, only: interpolateMT
      use source_location, only: locate_sources, rotate_mt
      use constants, only: IMAIN, DEBUG, &
                          dmom, dlat, dlon, ddep, dhdur, dcmt, &
                          partialnames
      use utils, only: gradient

      ! Inout
      type(t_GF), intent(in) :: GF
      type(t_source), dimension(1), intent(inout) :: sources ! Sources are located
      double precision, dimension(:,:,:), allocatable, intent(out) :: synt
      double precision, dimension(:,:,:,:), allocatable, intent(out) :: dp
      integer, intent(in) :: itypsokern

    end subroutine

  end interface get_seismograms

  ! Write seismograms interface
  interface write_seismograms
    module subroutine write_seismograms(&
      GF_filename, source_filename, output_dir, &
      OUTPUT_SEISMOS_SAC_ALPHANUM, &
      OUTPUT_SEISMOS_SAC_BINARY)
      character(len=*), intent(in) :: GF_filename, source_filename, output_dir
      logical, intent(in) :: OUTPUT_SEISMOS_SAC_ALPHANUM
      logical, intent(in) :: OUTPUT_SEISMOS_SAC_BINARY
    end subroutine

    module subroutine write_seismograms_sdp(&
      GF_filename, source_filename, output_dir, &
      itypsokern, &
      OUTPUT_SEISMOS_SAC_ALPHANUM, &
      OUTPUT_SEISMOS_SAC_BINARY)
      character(len=*), intent(in) :: GF_filename, source_filename, output_dir
      integer, intent(in) :: itypsokern
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