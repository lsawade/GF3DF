module gf3d_write_seismograms
  implicit none
  private
  public :: write_seismograms

contains

  subroutine write_seismograms(&
    GF_filename, source_filename, output_dir, &
    OUTPUT_SEISMOS_SAC_ALPHANUM, &
    OUTPUT_SEISMOS_SAC_BINARY)

    use gf3d_get_seismograms, only: get_seismograms
    use setup_source_location, only: setup_point_search_arrays
    use sac, only: write_output_SAC
    use constants, only: NCHANNELS, orientation, channels, MAX_STRING_LEN, IMAIN
    use sources, only: read_cmt, t_source
    use gf, only: t_GF, read_GF

    ! In
    character(len=*), intent(in) :: GF_filename, source_filename, output_dir
    logical, intent(in) :: OUTPUT_SEISMOS_SAC_ALPHANUM
    logical, intent(in) :: OUTPUT_SEISMOS_SAC_BINARY

    ! Local
    type(t_GF) :: GF
    type(t_source), dimension(:), allocatable :: sources
    character(len=MAX_STRING_LEN) :: sisname
    character(len=MAX_STRING_LEN) :: model = "GLAD-M25"
    integer :: icomp, k
    double precision, dimension(:,:,:), allocatable :: seismograms


    ! Read Green Function file
    GF = read_GF(GF_filename)

    ! Setup KDTree
    call setup_point_search_arrays(GF)

    ! Read cmt solution
    sources = read_cmt(source_filename)

    ! Extract seismograms
    seismograms = get_seismograms(GF, sources)

    do k=1,size(GF%displacement,1)



      do icomp=1,NCHANNELS
        write(IMAIN,*)  "Writing ", trim(GF%networks(k)), ".", trim(GF%stations(k)), ".", trim(channels(icomp))
        write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
        trim(GF%networks(k)), trim(GF%stations(k)), trim(channels(icomp))

        call write_output_SAC(&
          seismograms(k,icomp,:), &
          orientation(icomp), &
          sisname, &
          channels(icomp), &
          sources(1)%year, &
          sources(1)%jda, &
          sources(1)%hour, &
          sources(1)%minute, &
          sources(1)%second, &
          sources(1)%time_shift, &
          dble(GF%tc), &
          sources(1)%eventname, &
          sources(1)%latitude, &
          sources(1)%longitude, &
          sources(1)%depth, &
          sources(1)%hdur, &
          GF%stations(k), &
          GF%networks(k), &
          GF%latitudes(k), &
          GF%longitudes(k), &
          dble(0.0), &
          GF%burials(k), &
          GF%dt, &
          dble(0.0), &
          GF%nsteps, &
          OUTPUT_SEISMOS_SAC_ALPHANUM, &
          OUTPUT_SEISMOS_SAC_BINARY, &
          model, &
          output_dir)
      enddo
    enddo

  end subroutine write_seismograms

end module gf3d_write_seismograms