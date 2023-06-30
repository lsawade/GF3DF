submodule (gf3d) seismograms_windows
! This submodule is modeled after the the function call that Göran sent me a
! while back I have to of course change some things due to the way it is
! implemented, but in all in all it should work well with the CMT software.
!
! call getcmtpd_K52_swrays(stlat,stlon,stdip,stazi,tfir,spsam, &
!                          it1smp,it2smp,pminswrays,pmaxswrays, &
!                          isrc,synt,itypsokern,pdarray(1,iptr), &
!                          mxptsfit)
!
! INPUT:
!   stlat .... = station latitude
!   stlon .... = station longitude
!   stdip .... = dip of sensor axis
!   stazi .... = azimuth of sensor axis
!   tfir  .... = start time of first sample w.r.t. the reference time (event time)
!   spsam .... = seconds per sample (1/sample rate)
!   it1smp ... = index of first sample wanted
!   it2smp ... = index of last sample wanted
!   pminswrays = minimum period wanted
!   pmaxswrays = maximum period wanted
!   isrc ..... = source index (if more than one earthquake included in synthetics)
!   itypsokern = code for type of kernels desired (see below)
!   mxptsfit . = number of points in seismogram array

! OUTPUT:
!   synt ..... = synthetic seismogram synt(mxptsfit) from previous iteration
!   pdarray .. = array of partial derivatives pdarray(mxptsfit,index)

! Description:
!   if itypsokern=0, the subroutine only returns the synthetic seismogram in synt
!      itypsokern=1, the subroutine returns the synthetic seismogram in synt, and
!   the 6 moment-tensor kernels in pdarray
!      itypsokern=2, the subroutine returns the synthetic seismogram in synt, the
!   6 moment-tensor kernels pdarray (index 1-6) and 4 centroid kernels in
!      pdarray (index 7-10)

  implicit none
contains

  module subroutine get_seismograms_winsta(&
    GF, sources, &
    istat, &
    t0, dt, it0, itf,  &
    superseismograms)
    ! Retrieves a set of three seismograms given the Green Function class
    ! a single source or set of sources, timing, and station information
    !
    ! INPUT:
    !   net ...... = network name
    !   sta ...... = station name
    !   t0 ....... = start time of first sample w.r.t. the reference time (event time)
    !   dt ....... = seconds per sample (1/sample rate)
    !   it0 ...... = index of first sample wanted
    !   itf ...... = index of last sample wanted
    !   itypsokern = code for type of kernels desired (see below)
    !
    ! OUTPUT:
    !   superseismograms = OUTPUT seismogram array
    use sources, only: t_source
    use gf, only: t_GF
    use stf, only: get_stf, stf_convolution, correct_hdur
    use interpolation, only: interpolateMT, spline1d_seismograms
    use source_location, only: locate_sources
    use constants, only: IMAIN, DEBUG

    ! In
    type(t_GF), intent(in) :: GF
    type(t_source), dimension(:), intent(inout) :: sources
    integer, intent(in) :: istat
    double precision, intent(in) :: t0
    double precision, intent(in) :: dt
    integer, intent(in) :: it0
    integer, intent(in) :: itf

    ! Local
    integer :: isource, i, j, k, iglob, istat
    double precision :: t0_stf, tc_stf, hdur_diff
    double precision, dimension(:), allocatable :: t, stf
    double precision, dimension(:), allocatable :: tq
    double precision, dimension(:,:,:), allocatable :: seismograms
    double precision, dimension(:,:,:), allocatable :: convolution
    double precision, dimension(:,:,:,:,:,:,:), allocatable :: displacement
    double precision, dimension(:,:,:), allocatable :: interseismograms
    ! Out
    double precision, dimension(:,:,:), allocatable, intent(out) :: superseismograms

    ! Allocate all arrays
    allocate(t(GF%nsteps),stf(GF%nsteps))
    allocate(tq(itf-it0))
    allocate(seismograms(1, 3, GF%nsteps))
    allocate(convolution(1, 3, GF%nsteps))
    allocate(displacement(1, 3, 3, GF%ngllx,GF%nglly,GF%ngllz,GF%nsteps))
    allocate(interseismograms(1, 3, itf-it0))
    allocate(superseismograms(1, 3, itf-it0))

    ! Locate sources
    ! Put this somewhere outside this function otherwise we have to locate
    ! the source for every station
    call locate_sources(GF, sources)

    ! Setup basic parameter for STF
    t0_stf = 0.d0
    tc_stf = 200.d0
    t(:) = t0_stf + ((/(i, i=1, GF%nsteps, 1)/)-1) * GF%dt

    ! Create vector for interpolation in time that starts at t0 with respect to
    ! tc of the Green function database
    tq(:) = GF%tc + t0 + ((/(i, i=it0, itf, 1)/)-1) * dt

    ! Initialize
    superseismograms(:,:,:) = 0.d0

    do isource=1,size(sources)

      ! Rezero
      seismograms(:,:,:) = 0.d0
      convolution(:,:,:) = 0.d0
      interseismograms(:,:,:) = 0.d0
      displacement(:,:,:,:,:,:,:) = 0.d0
      stf(:) = 0.d0

      ! Get STF
      hdur_diff = correct_hdur(sources(isource)%hdur, GF%hdur)
      call get_stf(t, tc_stf, hdur_diff, int(GF%nsteps, kind=4), stf)

      ! Get displacement array
      do k=1,GF%ngllz
        do j=1,GF%nglly
          do i=1,GF%ngllx
            iglob = GF%ibool(i,j,k, sources(isource)%ispec)
            displacement(1,:,:,i,j,k, :) = GF%displacement(istat,:,:,iglob,:)
          enddo
        enddo
      enddo

      ! Interpolate displacement array to seismograms
      call interpolateMT(&
        displacement, 1, &
        GF%nsteps, GF%ngllx, GF%nglly, GF%ngllz, GF%xigll, GF%yigll, GF%zigll, &
        sources(isource)%Mxx, sources(isource)%Myy, sources(isource)%Mzz, &
        sources(isource)%Mxy, sources(isource)%Mxz, sources(isource)%Myz, &
        sources(isource)%xi, sources(isource)%eta,sources(isource)%gamma, &
        sources(isource)%xix, sources(isource)%xiy, sources(isource)%xiz, &
        sources(isource)%etax, sources(isource)%etay, sources(isource)%etaz, &
        sources(isource)%gammax, sources(isource)%gammay, sources(isource)%gammaz, &
        seismograms)

      ! Convolve seismogram array with STF
      convolution(:,:,:) = stf_convolution(seismograms, stf, GF%dt, tc_stf - sources(isource)%time_shift)

      ! Note that I'm reassigning seismograms here
      call spline1d_seismograms(t, convolution, tq, interseismograms)

      ! ADding to super seismogram array
      superseismograms(:,:,:) = superseismograms(:,:,:) + interseismograms(:,:,:)

      if (DEBUG) write(IMAIN,*) 'Max disp: ', maxval(displacement)
      if (DEBUG) write(IMAIN,*) 'Max seis: ', maxval(seismograms)
      if (DEBUG) write(IMAIN,*) 'Max stf:  ', maxval(stf)
      if (DEBUG) write(IMAIN,*) 'Max conv: ', maxval(convolution)
      if (DEBUG) write(IMAIN,*) 'Max supe: ', maxval(superseismograms)
    enddo


  end subroutine get_seismograms_winsta


  module subroutine write_seismograms_winsta(&
    GF_filename, source_filename, output_dir, &
    network, station, &
    t0, dt, it0, itf,  &
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
    character(len=*), intent(in) :: network, station
    logical, intent(in) :: OUTPUT_SEISMOS_SAC_ALPHANUM
    logical, intent(in) :: OUTPUT_SEISMOS_SAC_BINARY
    double precision, intent(in) :: t0
    double precision, intent(in) :: dt
    integer, intent(in) :: it0
    integer, intent(in) :: itf

    ! Local
    type(t_GF) :: GF
    type(t_source), dimension(:), allocatable :: sources
    character(len=MAX_STRING_LEN) :: sisname
    character(len=MAX_STRING_LEN) :: model = "GLAD-M25"
    integer :: i, icomp, k, istat
    double precision, dimension(:,:,:), allocatable :: seismograms


    ! Read Green Function file
    GF = read_GF(GF_filename)

    ! Setup KDTree
    call setup_point_search_arrays(GF)

    ! Read cmt solution
    sources = read_cmt(source_filename)

    ! Find station
    do i=1, size(GF%stations)

      if ((trim(network) == GF%networks(i)) .and. (trim(station) == GF%stations(i))) then
        istat = i
        exit
      else if (i==size(GF%stations)) then
        write (*,*) "Did not find ", trim(network), ".", trim(station), "in subset."
        stop "Went through all stations but requested one not found."
      endif

    end do

    ! Extract seismograms
    call get_seismograms_winsta(GF, sources, istat, t0, dt, it0, itf, seismograms)

    ! Write the three channels
    do icomp=1,NCHANNELS
       write(IMAIN,*)  "Writing ", trim(GF%networks(istat)), ".", trim(GF%stations(istat)), ".", trim(channels(icomp))
        write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
        trim(GF%networks(istat)), trim(GF%stations(istat)), trim(channels(icomp))

      call write_output_SAC(&
        seismograms(1,icomp,:), &
        orientation(icomp), &
        sisname, &
        channels(icomp), &
        sources(1)%year, &
        sources(1)%jda, &
        sources(1)%hour, &
        sources(1)%minute, &
        sources(1)%second, &
        sources(1)%time_shift, &
        - t0 - dt*it0, &
        sources(1)%eventname, &
        sources(1)%latitude, &
        sources(1)%longitude, &
        sources(1)%depth, &
        sources(1)%hdur, &
        station, &
        network, &
        GF%latitudes(istat), &
        GF%longitudes(istat), &
        dble(0.0), &
        GF%burials(istat), &
        dt, &
        dble(0.0), &
        int(itf-it0, kind=8), &
        OUTPUT_SEISMOS_SAC_ALPHANUM, &
        OUTPUT_SEISMOS_SAC_BINARY, &
        model, &
        output_dir)
    enddo


  end subroutine write_seismograms_winsta


end submodule seismograms_windows