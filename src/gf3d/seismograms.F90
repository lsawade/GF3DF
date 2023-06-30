submodule (gf3d) seismograms

  implicit none


  interface interpolate_source
    module subroutine interpolate_source(GF, source, seismograms)
      use gf, only: t_GF
      use sources, only: t_source
      use interpolation, only: interpolateMT
      type(t_GF), intent(in) :: GF
      type(t_source), intent(in) :: source
      double precision, dimension(:,:,:) :: seismograms
    end subroutine interpolate_source
  end interface interpolate_source

contains

  module subroutine get_sdp(GF, sources, synt, dp, itypsokern)

    use gf, only: t_GF
    use sources, only: t_source
    use stf, only: get_stf, stf_convolution, correct_hdur
    use interpolation, only: interpolateMT
    use source_location, only: locate_sources, rotate_mt
    use constants, only: IMAIN, DEBUG, &
                         dmom, dlat, dlon, ddep, dhdur, dcmt

    ! In
    type(t_GF), intent(in) :: GF
    type(t_source), dimension, intent(in) :: sources
    integer :: itypsokern

    ! Local
    integer :: i, j, k, iglob, Ndp
    double precision :: t0, tc, hdur_diff
    double precision, dimension(6) :: tM
    double precision, dimension(6,6) :: M
    type(t_source), dimension(1) :: dlat_source_p, dlat_source_m
    type(t_source), dimension(1) :: dlon_source_p, dlon_source_m
    type(t_source), dimension(1) :: dlon_source_p, dlon_source_m
    type(t_source), dimension(1) :: dhdur_source_p, dhdur_source_m
    double precision :: Mxx, Myy, Mzz, Mxy, Mxz, Myz
    double precision, dimension(:), allocatable :: t, stf, stf_p, stf_m
    double precision, dimension(:,:,:), allocatable :: convolution
    double precision, dimension(:,:,:,:,:,:,:), allocatable :: displacement
    double precision, dimension(:,:,:), allocatable :: seismograms

    ! Out
    double precision, dimension(:,:,:), allocatable, intent(out) :: synt
    double precision, dimension(:,:,:,:), allocatable, intent(out) :: dp

    ! Check which one of the partials to compute
    if (itypsokern == 1) then
      Ndp = 6
    elseif (itypsokern == 2) then
      Ndp = 10
    elseif (itypsokern == 3) then
      Ndp = 11
    else
      write(*,*) "Number of partials not implemented."
    endif

    ! Allocate all arrays
    allocate(t(GF%nsteps),stf(GF%nsteps))
    allocate(seismograms(size(GF%displacement,1), 3, GF%nsteps))
    allocate(convolution(size(GF%displacement,1), 3, GF%nsteps))
    allocate(displacement(size(GF%displacement,1), 3, 3, GF%ngllx,GF%nglly,GF%ngllz,GF%nsteps))
    allocate(synt(size(GF%displacement,1), 3, GF%nsteps))

    ! Allocate the array for partials
    allocate(dp(size(GF%displacement,1), 3, GF%nsteps))

    if (size(sources) .ne. 1) then
      stop "This function expects a single source. Check your sources array"
    endif

    ! Create base STF
    t0 = 0.d0
    tc = 200.d0
    t(:) = t0 + ((/(i, i=1, GF%nsteps, 1)/)-1) * GF%dt

    ! Get STF
    hdur_diff = correct_hdur(sources(isource)%hdur, GF%hdur)
    call get_stf(t, tc, hdur_diff, int(GF%nsteps, kind=4), stf)

    ! Get synthetics
    ! Rezero
    seismograms(:,:,:) = 0.d0
    convolution(:,:,:) = 0.d0
    displacement(:,:,:,:,:,:,:) = 0.d0

    ! Locate sources
    call locate_sources(GF, sources)

    ! Get synthetics

    ! Convolve seismogram array with STF
    convolution(:,:,:) = stf_convolution(seismograms, stf, GF%dt, tc - sources(isource)%time_shift)

    synt(:,:,:) =  convolution(:,:,:)

    if (DEBUG) write(IMAIN,*) 'Max disp: ', maxval(displacement)
    if (DEBUG) write(IMAIN,*) 'Max seis: ', maxval(seismograms)
    if (DEBUG) write(IMAIN,*) 'Max stf: ', maxval(stf)
    if (DEBUG) write(IMAIN,*) 'Max conv: ', maxval(convolution)
    if (DEBUG) write(IMAIN,*) 'Max supe: ', maxval(superseismograms)


    if (itypsokern > 0) then

      do i=1,6
        ! Copy source
        tM(:) = 0.d0
        tM(i) = dmom

        ! set all MT components to 0
        ! Assign moment tensor to array for simpler
        call rotate_mt(%
          tempsource%latitude, tempsource%longitude, &
          tM(1), tM(2), tM(3), tM(4), tM(5), tM(6), &
          M(i,1), M(i,2), M(i,3), M(i,4), M(i,5), M(i,6))
      enddo

    endif

    ! Create forward and back second order FD sources
    if (itypsokern > 1) then

      ! Latitude perturbation
      dlat_source_p(:) = sources(:)
      dlat_source_m(:) = sources(:)

      ! Longitude perturbation
      dlon_source_p(:) = sources(:)
      dlon_source_m(:) = sources(:)

      ! Depth perturbation
      ddep_source_p(:) = sources(:)
      ddep_source_m(:) = sources(:)

    endif

    !
    if (itypsokern > 2) then

      ! Half duration  perturbation
      dhdur_source_p(:) = sources(:)
      dhdur_source_m(:) = sources(:)
      allocate(stf_p(GF%nsteps), stf_m(GF%nsteps))

      ! Get STF
      hdur_diff = correct_hdur(sources(isource)%hdur, GF%hdur)
      call get_stf(t, tc, hdur_diff+dhdur, int(GF%nsteps, kind=4), stf_p)
      call get_stf(t, tc, hdur_diff-dhdur, int(GF%nsteps, kind=4), stf_m)

    endif



  end subroutine get_sdp


  module subroutine get_seismograms(GF, sources, superseismograms)

    use gf, only: t_GF
    use sources, only: t_source
    use stf, only: get_stf, stf_convolution, correct_hdur
    ! use interpolation, only: interpolateMT
    use source_location, only: locate_sources
    use constants, only: IMAIN, DEBUG

    ! In
    type(t_GF), intent(in) :: GF
    type(t_source), dimension(:), intent(inout) :: sources

    ! Local
    integer :: isource, i, j, k, iglob
    double precision :: t0, tc, hdur_diff
    double precision, dimension(:), allocatable :: t, stf
    double precision, dimension(:,:,:), allocatable :: seismograms
    double precision, dimension(:,:,:), allocatable :: convolution
    double precision, dimension(:,:,:,:,:,:,:), allocatable :: displacement

    ! Out
    double precision, dimension(:,:,:), allocatable, intent(out) :: superseismograms

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
    t(:) = t0 + ((/(i, i=1, GF%nsteps, 1)/)-1) * GF%dt

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

      ! Uses the sources located element coordinates and jacobians to
      ! Interpolate the traces
      call interpolate_source(GF, sources(isource), seismograms)

      ! Convolve seismogram array with STF
      convolution(:,:,:) = stf_convolution(seismograms, stf, GF%dt, tc - sources(isource)%time_shift)

      superseismograms(:,:,:) = superseismograms(:,:,:) + convolution(:,:,:)

      if (DEBUG) write(IMAIN,*) 'Max disp: ', maxval(displacement)
      if (DEBUG) write(IMAIN,*) 'Max seis: ', maxval(seismograms)
      if (DEBUG) write(IMAIN,*) 'Max stf: ', maxval(stf)
      if (DEBUG) write(IMAIN,*) 'Max conv: ', maxval(convolution)
      if (DEBUG) write(IMAIN,*) 'Max supe: ', maxval(superseismograms)
    enddo


  end subroutine get_seismograms



  module subroutine write_seismograms(&
    GF_filename, source_filename, output_dir, &
    OUTPUT_SEISMOS_SAC_ALPHANUM, &
    OUTPUT_SEISMOS_SAC_BINARY)

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
    call get_seismograms(GF, sources, seismograms)

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
          sources(1)%time_shift, & ! It's important to note that this has 0 effect!
          dble(GF%tc), & ! + sources(1)%time_shift , & ! Here we add the source time shift!
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



  ! module subroutine get_all_seismograms(&
  !   GF_filename, source_filename, output_dir, &
  !   OUTPUT_SEISMOS_SAC_ALPHANUM, &
  !   OUTPUT_SEISMOS_SAC_BINARY)

  !   use setup_source_location, only: setup_point_search_arrays
  !   use sac, only: write_output_SAC
  !   use constants, only: NCHANNELS, orientation, channels, MAX_STRING_LEN, IMAIN
  !   use sources, only: read_cmt, t_source
  !   use gf, only: t_GF, read_GF

  !   ! In
  !   character(len=*), intent(in) :: GF_filename, source_filename, output_dir
  !   logical, intent(in) :: OUTPUT_SEISMOS_SAC_ALPHANUM
  !   logical, intent(in) :: OUTPUT_SEISMOS_SAC_BINARY

  !   ! Local
  !   type(t_GF) :: GF
  !   type(t_source), dimension(:), allocatable :: sources
  !   character(len=MAX_STRING_LEN) :: sisname
  !   character(len=MAX_STRING_LEN) :: model = "GLAD-M25"
  !   integer :: icomp, k
  !   double precision, dimension(:,:,:), allocatable :: seismograms


  !   ! Read Green Function file
  !   GF = read_GF(GF_filename)

  !   ! Setup KDTree
  !   call setup_point_search_arrays(GF)

  !   ! Read cmt solution
  !   sources = read_cmt(source_filename)

  !   ! Extract seismograms
  !   call get_seismograms(GF, sources, seismograms)

  !   do k=1,size(GF%displacement,1)



  !     do icomp=1,NCHANNELS
  !       write(IMAIN,*)  "Writing ", trim(GF%networks(k)), ".", trim(GF%stations(k)), ".", trim(channels(icomp))
  !       write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
  !       trim(GF%networks(k)), trim(GF%stations(k)), trim(channels(icomp))

  !       call write_output_SAC(&
  !         seismograms(k,icomp,:), &
  !         orientation(icomp), &
  !         sisname, &
  !         channels(icomp), &
  !         sources(1)%year, &
  !         sources(1)%jda, &
  !         sources(1)%hour, &
  !         sources(1)%minute, &
  !         sources(1)%second, &
  !         sources(1)%time_shift, & ! It's important to note that this has 0 effect!
  !         dble(GF%tc), & ! + sources(1)%time_shift , & ! Here we add the source time shift!
  !         sources(1)%eventname, &
  !         sources(1)%latitude, &
  !         sources(1)%longitude, &
  !         sources(1)%depth, &
  !         sources(1)%hdur, &
  !         GF%stations(k), &
  !         GF%networks(k), &
  !         GF%latitudes(k), &
  !         GF%longitudes(k), &
  !         dble(0.0), &
  !         GF%burials(k), &
  !         GF%dt, &
  !         dble(0.0), &
  !         GF%nsteps, &
  !         OUTPUT_SEISMOS_SAC_ALPHANUM, &
  !         OUTPUT_SEISMOS_SAC_BINARY, &
  !         model, &
  !         output_dir)
  !     enddo
  !   enddo

  ! end subroutine get_all_seismograms


end submodule seismograms