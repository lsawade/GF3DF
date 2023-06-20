
module gf3d

  use ctypes, only: t_GF, t_source
  use io, only: read_GF, free_GF, read_cmt, print_GF, print_source, write_output_SAC
  use utils, only: throwerror, scaleM
  use setup_source_location, only: setup_point_search_arrays
  use source_location, only: locate_sources
  use interpolation, only: interpolateMT
  use hdf5

  private
  public :: &
    read_GF, read_cmt, &
    print_GF, print_source, &
    t_GF, &
    t_source, &
    free_GF, &
    throwerror, &
    setup_point_search_arrays, &
    locate_sources, &
    scaleM, &
    julian_day, &
    interpolateMT, &
    write_output_SAC, &
    write_seismograms


contains

  subroutine write_seismograms(&
    GF_filename, source_filename, output_dir, &
    OUTPUT_SEISMOS_SAC_ALPHANUM, &
    OUTPUT_SEISMOS_SAC_BINARY)

    use constants, only: NCHANNELS, orientation, channels, MAX_STRING_LEN
    use io, only: read_cmt, read_GF
    use ctypes, only: t_GF, t_source

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

    ! Read cmt solution
    sources = read_cmt(source_filename)

    ! Extract seismograms
    seismograms = get_seismograms(GF, sources)

    do k=1,size(GF%displacement,1)

      write(*,*) trim(GF%networks(k)), ".", trim(GF%stations(k))

      do icomp=1,NCHANNELS

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


  function get_seismograms(GF, sources) result(superseismograms)

    use ctypes, only: t_GF, t_source
    use stf, only: get_stf, stf_convolution, correct_hdur
    use interpolation, only: interpolateMT

    ! In
    type(t_GF), intent(in) :: GF
    type(t_source), dimension(:) :: sources

    ! Local
    integer :: isource
    double precision :: t0, tc, hdur_diff
    double precision, dimension(:), allocatable :: t, stf
    double precision, dimension(:,:,:), allocatable :: seismograms
    double precision, dimension(:,:,:), allocatable :: convolution
    double precision, dimension(:,:,:,:,:,:,:), allocatable :: displacement

    ! Out
    double precision, dimension(:,:,:), allocatable :: superseismograms

    ! Allocate all arrays
    allocate(t(GF%nsteps),stf(GF%nsteps))
    allocate(seismograms(size(GF%displacement,1), 3, GF%nsteps))
    allocate(convolution(size(GF%displacement,1), 3, GF%nsteps))
    allocate(superseismograms(size(GF%displacement,1), 3, GF%nsteps))
    allocate(displacement(size(GF%displacement,1), 3, 3, GF%ngllx,GF%nglly,GF%ngllz,GF%nsteps))

    ! Locate sources
    call locate_sources(GF, sources)

    ! Setup basic parameter for STF
    t0 = 0.d0
    tc = 200.d0
    t(:) = t0 + ((/(I, I=1, GF%nsteps, 1)/)-1) * GF%dt

    ! Initialize
    superseismograms(:,:,:) = 0.d0

    do isource=1,size(sources)

      ! Rezero
      seismograms(:,:,:) = 0.d0
      convolution(:,:,:) = 0.d0
      displacement(:,:,:,:,:,:,:) = 0.d0
      stf(:) = 0.d0

      ! Get STF
      hdur_diff = correct_hdur(sources(isource)%hdur, GF%hdur)
      call get_stf(t, tc, hdur_diff, int(GF%nsteps, kind=4), stf)

      ! Get displacement array
      do k=1,GF%ngllz
        do j=1,GF%nglly
          do i=1,GF%ngllx
            iglob = GF%ibool(i,j,k, sources(isource)%ispec)
            displacement(:,:,:,i,j,k, :) = GF%displacement(:,:,:,iglob,:)
          enddo
        enddo
      enddo

      ! Interpolate displacement array to seismograms
      call interpolateMT(&
        displacement, size(GF%displacement,1), &
        GF%nsteps, GF%ngllx, GF%nglly, GF%ngllz, GF%xigll, GF%yigll, GF%zigll, &
        sources(isource)%Mxx, sources(isource)%Myy, sources(isource)%Mzz, &
        sources(isource)%Mxy, sources(isource)%Mxz, sources(isource)%Myz, &
        sources(isource)%xi, sources(isource)%eta,sources(isource)%gamma, &
        sources(isource)%xix, sources(isource)%xiy, sources(isource)%xiz, &
        sources(isource)%etax, sources(isource)%etay, sources(isource)%etaz, &
        sources(isource)%gammax, sources(isource)%gammay, sources(isource)%gammaz, &
        seismograms)

      ! Convolve seismogram array with STF
      convolution(:,:,:) = stf_convolution(seismograms, stf, GF%dt, tc)

      superseismograms(:,:,:) = superseismograms(:,:,:) + convolution(:,:,:)

      write (*,*) 'Max disp: ', maxval(displacement)
      write (*,*) 'Max seis: ', maxval(seismograms)
      write (*,*) 'Max stf: ', maxval(stf)
      write (*,*) 'Max conv: ', maxval(convolution)
      write (*,*) 'Max supe: ', maxval(superseismograms)
    enddo


  end function get_seismograms


end module gf3d