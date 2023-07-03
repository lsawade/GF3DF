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
                         dmom, dlat, dlon, ddep, dhdur, dcmt, &
                         partialnames
    use utils, only: gradient

    ! In
    type(t_GF), intent(in) :: GF
    type(t_source), dimension(1), intent(inout) :: sources
    integer :: itypsokern

    ! Local
    integer :: i, j, k, iglob, Ndp
    double precision :: t0, tc, hdur_diff
    double precision, dimension(6) :: tM
    type(t_source), dimension(6) :: dmom_source
    type(t_source), dimension(1) :: dlat_source_p, dlat_source_m
    type(t_source), dimension(1) :: dlon_source_p, dlon_source_m
    type(t_source), dimension(1) :: ddep_source_p, ddep_source_m
    type(t_source), dimension(1) :: dhdur_source_p, dhdur_source_m
    double precision :: Mxx, Myy, Mzz, Mxy, Mxz, Myz
    double precision, dimension(:), allocatable :: t, stf, stf_p, stf_m
    double precision, dimension(:,:,:), allocatable :: convolution
    double precision, dimension(:,:,:,:,:,:,:), allocatable :: displacement
    double precision, dimension(:,:,:), allocatable :: seismograms
    double precision, dimension(:,:,:), allocatable :: seismograms2

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
    allocate(convolution(size(GF%displacement,1), 3, GF%nsteps))
    allocate(seismograms(size(GF%displacement,1), 3, GF%nsteps))
    allocate(seismograms2(size(GF%displacement,1), 3, GF%nsteps))
    allocate(synt(size(GF%displacement,1), 3, GF%nsteps))

    ! Allocate the array for partials
    allocate(dp(Ndp, size(GF%displacement,1), 3, GF%nsteps))

    if (DEBUG) write(IMAIN, *) "Allocated everything."

    if (size(sources) .ne. 1) then
      stop "This function expects a single source. Check your sources array"
    endif

    ! Create base STF
    t0 = 0.d0
    tc = 200.d0
    t(:) = t0 + ((/(i, i=1, GF%nsteps, 1)/)-1) * GF%dt

    ! Get STF
    hdur_diff = correct_hdur(sources(1)%hdur, GF%hdur)
    call get_stf(t, tc, hdur_diff, int(GF%nsteps, kind=4), stf)

    if (DEBUG) write(IMAIN, *) "Created bas STF."

    ! Get synthetics
    ! Rezero
    seismograms(:,:,:) = 0.d0
    convolution(:,:,:) = 0.d0
    synt(:,:,:) = 0.d0
    dp(:,:,:,:) = 0.d0

    ! Locate sources
    call locate_sources(GF, sources)

    if (DEBUG) write(IMAIN, *) "Located base source."

    ! Get synthetics
    call interpolate_source(GF, sources(1), seismograms)

    ! Convolve seismogram array with STF
    convolution(:,:,:) = stf_convolution(seismograms, stf, GF%dt, tc - sources(1)%time_shift)

    if (DEBUG) write(IMAIN, *) "Convolved base seismograms."

    synt(:,:,:) = convolution(:,:,:)

    if (DEBUG) write(IMAIN,*) 'Max seis: ', maxval(seismograms)
    if (DEBUG) write(IMAIN,*) 'Max stf:  ', maxval(stf)
    if (DEBUG) write(IMAIN,*) 'Max conv: ', maxval(convolution)
    if (DEBUG) write(IMAIN,*) 'Max supe: ', maxval(synt)

    if (itypsokern > 0) then

      do i=1,6
        ! Copy source
        tM(:) = 0.d0
        tM(i) = dmom

        dmom_source(i) = sources(1)

        ! set all MT components to 0
        ! Assign moment tensor to array for simpler
        call rotate_mt(&
          sources(1)%latitude, sources(1)%longitude, &
          tM(1), tM(2), tM(3), tM(4), tM(5), tM(6), &
          dmom_source(i)%Mxx, dmom_source(i)%Myy, dmom_source(i)%Mzz, &
          dmom_source(i)%Mxy, dmom_source(i)%Mxz, dmom_source(i)%Myz)

        ! Get synthetics
        call interpolate_source(GF, dmom_source(i), seismograms)

        ! Convolve seismogram array with STF
        convolution(:,:,:) = stf_convolution(seismograms, stf, GF%dt, tc - sources(1)%time_shift)

        dp(i,:,:,:) =  convolution(:,:,:)/dmom

        if (DEBUG) write(IMAIN, *) "Done with ", trim(partialnames(i))

      enddo

    endif

    ! Create forward and back second order FD sources
    if (itypsokern > 1) then

      ! ==================== Latitude

      ! Latitude perturbation
      dlat_source_p(1) = sources(1)
      dlat_source_m(1) = sources(1)

      ! Perturb latitude locations
      dlat_source_p(1)%latitude = dlat_source_p(1)%latitude + dlat
      dlat_source_m(1)%latitude = dlat_source_p(1)%latitude - dlat

      ! Locate new locations
      call locate_sources(GF, dlat_source_p)
      call locate_sources(GF, dlat_source_m)

      ! Get synthetics
      ! Interpolate displacement seismograms for both locations
      call interpolate_source(GF, dlat_source_p(1), seismograms)
      call interpolate_source(GF, dlat_source_m(1), seismograms2)

      ! We can convolove afterwards since convolution is distributive
      ! Convolve seismogram array with STF
      convolution(:,:,:) = stf_convolution(&
        (seismograms-seismograms2)/(2*dlat), stf, &
        GF%dt, tc - sources(1)%time_shift)

      dp(7,:,:,:) =  convolution(:,:,:)

      if (DEBUG) write(IMAIN, *) "Done with ", trim(partialnames(7))

      ! ==================== Longitude

      ! Longitude perturbation
      dlon_source_p(:) = sources(:)
      dlon_source_m(:) = sources(:)

      ! Perturb latitude locations
      dlon_source_p(1)%longitude = dlon_source_p(1)%longitude + dlon
      dlon_source_m(1)%longitude = dlon_source_p(1)%longitude - dlon

      ! Locate new locations
      call locate_sources(GF, dlon_source_p)
      call locate_sources(GF, dlon_source_m)

      ! Interpolate displacement seismograms for both locations
      call interpolate_source(GF, dlon_source_p(1), seismograms)
      call interpolate_source(GF, dlon_source_m(1), seismograms2)

      convolution(:,:,:) = stf_convolution(&
        (seismograms-seismograms2)/(2*dlon), stf, &
        GF%dt, tc - sources(1)%time_shift)

      dp(8,:,:,:) =  convolution(:,:,:)

      if (DEBUG) write(IMAIN, *) "Done with ", trim(partialnames(8))

      ! =================== Depth

      ! Depth perturbation
      ddep_source_p(:) = sources(:)
      ddep_source_m(:) = sources(:)

      ! Perturb latitude locations
      ddep_source_p(1)%longitude = ddep_source_p(1)%longitude + ddep
      ddep_source_m(1)%longitude = ddep_source_p(1)%longitude - ddep

      ! Locate new locations
      call locate_sources(GF, ddep_source_p)
      call locate_sources(GF, ddep_source_m)

      ! Interpolate displacement seismograms for both locations
      call interpolate_source(GF, ddep_source_p(1), seismograms)
      call interpolate_source(GF, ddep_source_m(1), seismograms2)

      convolution(:,:,:) = stf_convolution(&
        (seismograms-seismograms2)/(2*ddep), stf, &
        GF%dt, tc - sources(1)%time_shift)

      dp(9,:,:,:) =  convolution(:,:,:)

      if (DEBUG) write(IMAIN, *) "Done with ", trim(partialnames(9))

      ! =================== Time

      call gradient(synt, Gf%dt, dp(10,:,:,:))

      dp(10,:,:,:) = (-1.d0) * dp(10,:,:,:)

      if (DEBUG) write(IMAIN, *) "Done with ", trim(partialnames(10))

    endif

    !
    if (itypsokern > 2) then

      ! Half duration  perturbation
      dhdur_source_p(:) = sources(:)
      dhdur_source_m(:) = sources(:)
      allocate(stf_p(GF%nsteps), stf_m(GF%nsteps))

      ! Get STF
      hdur_diff = correct_hdur(sources(1)%hdur, GF%hdur)
      call get_stf(t, tc, hdur_diff+dhdur, int(GF%nsteps, kind=4), stf_p)
      call get_stf(t, tc, hdur_diff-dhdur, int(GF%nsteps, kind=4), stf_m)

      ! compute partial derivatives with respect to half duration
      dp(11,:,:,:) = (&
        stf_convolution(seismograms, stf_p, &
                        GF%dt, tc - sources(1)%time_shift) - &
        stf_convolution(seismograms, stf_m, &
                        GF%dt, tc - sources(1)%time_shift))/dhdur

      if (DEBUG) write(IMAIN, *) "Done with ", trim(partialnames(11))
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


  module subroutine write_seismograms_sdp(&
    GF_filename, source_filename, output_dir, &
    itypsokern, &
    OUTPUT_SEISMOS_SAC_ALPHANUM, &
    OUTPUT_SEISMOS_SAC_BINARY)

    use setup_source_location, only: setup_point_search_arrays
    use sac, only: write_output_SAC
    use constants, only: NCHANNELS, orientation, channels, MAX_STRING_LEN, IMAIN, partialnames
    use sources, only: read_cmt, t_source
    use gf, only: t_GF, read_GF

    ! In
    character(len=*), intent(in) :: GF_filename, source_filename, output_dir
    integer, intent(in) :: itypsokern
    logical, intent(in) :: OUTPUT_SEISMOS_SAC_ALPHANUM
    logical, intent(in) :: OUTPUT_SEISMOS_SAC_BINARY

    ! Local
    type(t_GF) :: GF
    type(t_source), dimension(:), allocatable :: sources
    character(len=MAX_STRING_LEN) :: sisname
    character(len=MAX_STRING_LEN) :: model = "GLAD-M25"
    integer :: icomp, k, ip, Ndp = 0
    double precision, dimension(:,:,:), allocatable :: synt
    double precision, dimension(:,:,:,:), allocatable :: dp

    ! Read Green Function file
    GF = read_GF(GF_filename)

    ! Setup KDTree
    call setup_point_search_arrays(GF)

    ! Read cmt solution
    sources = read_cmt(source_filename)

    ! Extract seismograms
    call get_sdp(GF, sources, synt, dp, itypsokern)

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


    do k=1,size(GF%displacement,1)

      do icomp=1,NCHANNELS
        write(IMAIN,*)  "Writing ", trim(GF%networks(k)), ".", trim(GF%stations(k)), ".", trim(channels(icomp))
        write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
        trim(GF%networks(k)), trim(GF%stations(k)), trim(channels(icomp))

        call write_output_SAC(&
          synt(k,icomp,:), &
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

        if (Ndp .ne. 0) then
          do ip=1,Ndp

            write(IMAIN,*)  "Writing ", trim(GF%networks(k)), ".", trim(GF%stations(k)), ".", &
                                        trim(channels(icomp)), ".", trim(partialnames(ip))
            write(sisname,"('/',a,'.',a,'.',a3,'.',a,'.sem')") &
            trim(GF%networks(k)), trim(GF%stations(k)), trim(channels(icomp)), trim(partialnames(ip))

          call write_output_SAC(&
            dp(ip, k,icomp,:), &
            orientation(icomp), &
            sisname, &
            channels(icomp), &
            sources(1)%year, &
            sources(1)%jda, &
            sources(1)%hour, &
            sources(1)%minute, &
            sources(1)%second, &
            sources(1)%time_shift, & ! This does not have an effect
            dble(GF%tc), & ! This has an effect
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
        endif

      enddo
    enddo

  end subroutine write_seismograms_sdp



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