!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!----
!----  locate_sources finds the correct position of the sources
!----  This code is heavily modified for cleanliness and simplicity
!----
module source_location

  private
  public :: locate_sources

contains

  subroutine locate_sources(GF, sources)

    use constants, only: R_PLANET, RHOAV, ZERO, &
                        SOURCE_DECAY_MIMIC_TRIANGLE, RADIANS_TO_DEGREES, &
                        R_UNIT_SPHERE, PI_OVER_TWO, PI, IMAIN, GRAV, &
                        DEGREES_TO_RADIANS, TWO_PI
    use topo, only: get_topo_bathy
    use spline_routines, only: spline_evaluation
    use utils, only: throwerror
    use rthetaphi_xyz, only: lat_2_geocentric_colat_dble, xyz_2_rthetaphi_dble, &
                             geocentric_2_geographic_dble
    use point_location, only: locate_point
    use sources, only: &
      t_source, &
      get_cmt_scalar_moment, &
      get_cmt_moment_magnitude, &
      get_cmt_moment_magnitude_from_M0

    ! use reduce
    use gf, only: t_GF

    implicit none

    type(t_GF), intent(in) :: GF
    type(t_source), dimension(:), intent(inout), target :: sources

    ! local parameters
    integer :: isource, NSOURCES

    ! Source parameter pointers
    double precision :: depth
    double precision :: lat
    double precision :: lon
    double precision, pointer :: Mrr
    double precision, pointer :: Mtt
    double precision, pointer :: Mpp
    double precision, pointer :: Mrt
    double precision, pointer :: Mrp
    double precision, pointer :: Mtp
    double precision, pointer :: Mxx
    double precision, pointer :: Myy
    double precision, pointer :: Mzz
    double precision, pointer :: Mxy
    double precision, pointer :: Mxz
    double precision, pointer :: Myz
    double precision, dimension(:,:), pointer :: nu

    double precision :: sint,cost,sinp,cosp

    double precision :: ell
    double precision :: elevation
    double precision :: r0,p20, radius

    double precision :: source_final_distance_max

    integer :: iorientation
    double precision :: stazi,stdip
    double precision :: n(3),thetan,phin

    double precision :: Mw, M0
    double precision :: f0,t0_ricker,scaleF,force_N
    double precision :: total_M0,total_Mw,total_force_N


    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating sources'
    write(IMAIN,*) '********************'
    write(IMAIN,*)



    total_M0 = 0.d0
    total_Mw = 0.d0
    total_force_N = 0.d0

    ! normalized source radius
    r0 = R_UNIT_SPHERE

    ! loop on all the sources
    ! gather source information in subsets to reduce memory requirements
    NSOURCES = size(sources)

    ! loop over subsets of sources
    do isource= 1, NSOURCES

      ! arrays to collect data

      ! initializes
      ! loop over sources within this subset

      ! source lat/lon in degrees
      lat = sources(isource)%latitude
      lon = sources(isource)%longitude

      ! limits longitude to [0.0,360.0]
      if (lon < 0.d0 ) lon = lon + 360.d0
      if (lon > 360.d0 ) lon = lon - 360.d0

      ! convert geographic latitude lat (degrees) to geocentric colatitude theta (radians)
      call lat_2_geocentric_colat_dble(lat,sources(isource)%theta)

      sources(isource)%phi = lon*DEGREES_TO_RADIANS

      call reduce(sources(isource)%theta,sources(isource)%phi)

      sint = sin(sources(isource)%theta)
      cost = cos(sources(isource)%theta)
      sinp = sin(sources(isource)%phi)
      cosp = cos(sources(isource)%phi)

      if (sources(isource)%force .eqv. .false.) then

        if (isource == 1) write(IMAIN, *) " Rotating moment tensor"

        ! get the moment tensor
        Mrr => sources(isource)%Mrr
        Mtt => sources(isource)%Mtt
        Mpp => sources(isource)%Mpp
        Mrt => sources(isource)%Mrt
        Mrp => sources(isource)%Mrp
        Mtp => sources(isource)%Mtp

        Mxx => sources(isource)%Mxx
        Myy => sources(isource)%Myy
        Mzz => sources(isource)%Mzz
        Mxy => sources(isource)%Mxy
        Mxz => sources(isource)%Mxz
        Myz => sources(isource)%Myz

        ! convert from a spherical to a Cartesian representation of the moment tensor
        Mxx = sint*sint*cosp*cosp*Mrr + cost*cost*cosp*cosp*Mtt + sinp*sinp*Mpp &
            + 2.0d0*sint*cost*cosp*cosp*Mrt - 2.0d0*sint*sinp*cosp*Mrp - 2.0d0*cost*sinp*cosp*Mtp

        Myy = sint*sint*sinp*sinp*Mrr + cost*cost*sinp*sinp*Mtt + cosp*cosp*Mpp &
            + 2.0d0*sint*cost*sinp*sinp*Mrt + 2.0d0*sint*sinp*cosp*Mrp + 2.0d0*cost*sinp*cosp*Mtp

        Mzz = cost*cost*Mrr + sint*sint*Mtt - 2.0d0*sint*cost*Mrt

        Mxy = sint*sint*sinp*cosp*Mrr + cost*cost*sinp*cosp*Mtt - sinp*cosp*Mpp &
            + 2.0d0*sint*cost*sinp*cosp*Mrt + sint*(cosp*cosp-sinp*sinp)*Mrp + cost*(cosp*cosp-sinp*sinp)*Mtp

        Mxz = sint*cost*cosp*Mrr - sint*cost*cosp*Mtt &
            + (cost*cost-sint*sint)*cosp*Mrt - cost*sinp*Mrp + sint*sinp*Mtp

        Myz = sint*cost*sinp*Mrr - sint*cost*sinp*Mtt &
            + (cost*cost-sint*sint)*sinp*Mrt + cost*cosp*Mrp - sint*cosp*Mtp


      endif

      ! record three components for each station
      do iorientation = 1,3

        nu => sources(isource)%nu

        !   North
        if (iorientation == 1) then
          stazi = 0.d0
          stdip = 0.d0
        !   East
        else if (iorientation == 2) then
          stazi = 90.d0
          stdip = 0.d0
        !   Vertical
        else if (iorientation == 3) then
          stazi = 0.d0
          stdip = - 90.d0
        else
          call throwerror(1, "Wrong orientation.")
        endif

        !   get the orientation of the seismometer
        thetan = (90.0d0+stdip)*DEGREES_TO_RADIANS
        phin = stazi*DEGREES_TO_RADIANS

        ! we use the same convention as in Harvard normal modes for the orientation

        !   vertical component
        n(1) = cos(thetan)
        !   N-S component
        n(2) = - sin(thetan)*cos(phin)
        !   E-W component
        n(3) = sin(thetan)*sin(phin)

        !   get the Cartesian components of n in the model: nu
        nu(iorientation,1) = n(1)*sint*cosp + n(2)*cost*cosp - n(3)*sinp
        nu(iorientation,2) = n(1)*sint*sinp + n(2)*cost*sinp + n(3)*cosp
        nu(iorientation,3) = n(1)*cost - n(2)*sint

      enddo


      ! point depth (in m)
      depth = sources(isource)%depth * 1000.0d0

      write(*,*) "MT depth"
      ! normalized source radius
      r0 = R_UNIT_SPHERE

      ! finds elevation of position
      if (GF%topography == 1) then
        call get_topo_bathy(lat,lon,elevation, GF%nx_bathy, GF%ny_bathy, GF%bathy, GF%resolution_topo_file)
        r0 = r0 + elevation/R_PLANET
      else if (GF%topography == 0) then
        write(IMAIN, *) "No topography applied."
      else
        call throwerror(1, "Wrong topography value.")
      endif

      ! ellipticity
      if (GF%ellipticity == 1) then
        ! this is the Legendre polynomial of degree two, P2(cos(theta)),
        ! see the discussion above eq (14.4) in Dahlen and Tromp (1998)
        p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)

        ! todo: check if we need radius or r0 for evaluation below...
        !       (receiver location routine takes r0)
        radius = r0 - depth/R_PLANET

        ! get ellipticity using spline evaluation
        call spline_evaluation(GF%rspl,GF%ellipticity_spline,GF%ellipticity_spline2, size(GF%ellipticity_spline),radius,ell)

        ! this is eq (14.4) in Dahlen and Tromp (1998)
        r0 = r0*(1.0d0-(2.0d0/3.0d0)*ell*p20)
      endif

      ! stores surface radius for info output
      sources(isource)%r0 = r0

      ! subtracts source depth (given in m)
      sources(isource)%r_target = sources(isource)%r0 - depth/R_PLANET

      ! compute the Cartesian position of the source
      sources(isource)%x_target = sources(isource)%r_target*sint*cosp
      sources(isource)%y_target = sources(isource)%r_target*sint*sinp
      sources(isource)%z_target = sources(isource)%r_target*cost

      write (*,*) "Anchors"
      write (*,*) GF%anchor_iax
      write (*,*) GF%anchor_iay
      write (*,*) GF%anchor_iaz
      ! locates best element and xi/eta/gamma interpolation values
      call locate_point(&
        sources(isource)%x_target, &
        sources(isource)%y_target, &
        sources(isource)%z_target, &
        sources(isource)%latitude, &
        sources(isource)%longitude, &
        sources(isource)%ispec, &
        GF%nspec, GF%ngllx, GF%nglly, GF%ngllz, GF%midx, GF%midy, GF%midz, &
        GF%ibool, GF%xyz(:,1),GF%xyz(:,2), GF%xyz(:,3), GF%xadj, GF%adjacency, &
        GF%xigll, GF%yigll, GF%zigll, GF%anchor_iax, GF%anchor_iay, GF%anchor_iaz, &
        .true., &
        sources(isource)%xi, &
        sources(isource)%eta, &
        sources(isource)%gamma, &
        sources(isource)%xix,sources(isource)%xiy,sources(isource)%xiz, &
        sources(isource)%etax,sources(isource)%etay,sources(isource)%etaz, &
        sources(isource)%gammax,sources(isource)%gammay,sources(isource)%gammaz, &
        sources(isource)%x, &
        sources(isource)%y, &
        sources(isource)%z, &
        sources(isource)%final_distance)


        ! source info
      write(IMAIN,*)
      write(IMAIN,*) 'source # ',isource
      write(IMAIN,*)
      write(IMAIN,*) '  source located in element ',sources(isource)%ispec
      write(IMAIN,*)
      ! different output for force point sources
      if (sources(isource)%force) then
        write(IMAIN,*) '  using force point source: '
        write(IMAIN,*) '    xi coordinate of source in that element: ',sources(isource)%xi
        write(IMAIN,*) '    eta coordinate of source in that element: ',sources(isource)%eta
        write(IMAIN,*) '    gamma coordinate of source in that element: ',sources(isource)%gamma

        write(IMAIN,*)
        write(IMAIN,*) '    component of direction vector in East direction: ', &
            sources(isource)%comp_dir_vect_source_E
        write(IMAIN,*) '    component of direction vector in North direction: ', &
            sources(isource)%comp_dir_vect_source_N
        write(IMAIN,*) '    component of direction vector in Vertical direction: ', &
            sources(isource)%comp_dir_vect_source_Z_UP

        !write(IMAIN,*) '  i index of source in that element: ',nint(xi_source(isource))
        !write(IMAIN,*) '  j index of source in that element: ',nint(eta_source(isource))
        !write(IMAIN,*) '  k index of source in that element: ',nint(gamma_source(isource))
        !write(IMAIN,*)
        !write(IMAIN,*) '  component direction: ',COMPONENT_FORCE_SOURCE
        write(IMAIN,*)
        write(IMAIN,*) '    nu1 = ',nu(1,:),'North'
        write(IMAIN,*) '    nu2 = ',nu(2,:),'East'
        write(IMAIN,*) '    nu3 = ',nu(3,:),'Vertical'
        write(IMAIN,*)
        write(IMAIN,*) '    at (x,y,z) coordinates = ', sources(isource)%x, &
          sources(isource)%y, sources(isource)%z
      else
        ! moment tensor
        write(IMAIN,*) '  using moment tensor source: '
        write(IMAIN,*) '    xi coordinate of source in that element: ',sources(isource)%xi
        write(IMAIN,*) '    eta coordinate of source in that element: ',sources(isource)%eta
        write(IMAIN,*) '    gamma coordinate of source in that element: ',sources(isource)%gamma

        write(IMAIN,*) '    at (x,y,z) coordinates = ', sources(isource)%x, &
          sources(isource)%y, sources(isource)%z
      endif
      write(IMAIN,*)

      ! source time function info
      write(IMAIN,*) '  source time function:'
      ! frequency/half-duration
      if (sources(isource)%force) then
        ! single point force
        ! prints frequency content for point forces
        select case(sources(isource)%force_stf)
        case (0)
          ! Gaussian
          write(IMAIN,*) '    using Gaussian source time function'
          write(IMAIN,*) '             half duration: ',sources(isource)%hdur,' seconds'
          write(IMAIN,*) '    Gaussian half duration: ',sources(isource)%hdur/SOURCE_DECAY_MIMIC_TRIANGLE,' seconds'
        case (1)
          ! Ricker
          write(IMAIN,*) '    using Ricker source time function'
          ! prints frequency content for point forces
          f0 = sources(isource)%hdur
          t0_ricker = 1.2d0/f0
          write(IMAIN,*)
          write(IMAIN,*) '    using a source of dominant frequency ',f0
          write(IMAIN,*) '    t0_ricker = ',t0_ricker,'tshift_src = ',sources(isource)%time_shift
          write(IMAIN,*)
          write(IMAIN,*) '    lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
          write(IMAIN,*) '    lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
          write(IMAIN,*)
          write(IMAIN,*) '    half duration in frequency: ',sources(isource)%hdur,' seconds**(-1)'
        case (2)
          ! Heaviside
          write(IMAIN,*) '    using (quasi) Heaviside source time function'
          write(IMAIN,*) '             half duration: ',sources(isource)%hdur,' seconds'
        case (3)
          ! Monochromatic
          write(IMAIN,*) '    using monochromatic source time function'
          ! prints frequency content for point forces
          f0 = sources(isource)%hdur
          write(IMAIN,*)
          write(IMAIN,*) '    using a source of period ',f0
          write(IMAIN,*)
          write(IMAIN,*) '    half duration in period: ',sources(isource)%hdur,' seconds'
        case (4)
          ! Gaussian by Meschede et al. (2011)
          write(IMAIN,*) '    using Gaussian source time function by Meschede et al. (2011), eq.(2)'
          write(IMAIN,*) '             tau: ',sources(isource)%hdur,' seconds'
        case default
          stop 'unsupported force_stf value!'
        end select

      else
        ! moment tensor
        write(IMAIN,*) '    using (quasi) Heaviside source time function'
        ! add message if source is a Heaviside
        write(IMAIN,*)
        write(IMAIN,*) '    half duration: ',sources(isource)%hdur,' seconds'
      endif
      write(IMAIN,*) '    time shift: ',sources(isource)%time_shift,' seconds'
      write(IMAIN,*)

      ! magnitude
      write(IMAIN,*) '  magnitude of the source:'
      if (sources(isource)%force) then
        ! single point force
        ! scale and non-dimensionalize the factor_force_source
        ! factor_force_source in FORCESOLUTION file is in Newton
        ! 1 Newton is 1 kg * 1 m / (1 second)^2
        scaleF = RHOAV * (R_PLANET**4) * PI*GRAV*RHOAV
        ! force in Newton
        force_N = sources(isource)%factor_force_source * scaleF
        ! adds to total force applied (sum over all force point sources)
        total_force_N = total_force_N + force_N

        write(IMAIN,*) '    force = ', sngl(force_N),'(Newton)' ! dimensionalized
      else
        ! moment-tensor
        sources(isource)%M0 = get_cmt_scalar_moment(Mxx,Myy,Mzz,Mxy,Mxz,Myz)
        sources(isource)%Mw =  get_cmt_moment_magnitude(Mxx,Myy,Mzz,Mxy,Mxz,Myz)
        ! adds to total moment
        total_M0 = total_M0 + M0
        total_Mw = get_cmt_moment_magnitude_from_M0(total_M0)

        write(IMAIN,*) '       scalar moment M0 = ', sources(isource)%M0,' dyne-cm'
        write(IMAIN,*) '    moment magnitude Mw = ', sources(isource)%Mw
      endif
      write(IMAIN,*)

      ! get latitude, longitude and depth of the source that will be used
      call xyz_2_rthetaphi_dble(sources(isource)%x, &
                                sources(isource)%y, &
                                sources(isource)%z, &
                                sources(isource)%r_found, &
                                sources(isource)%theta, &
                                sources(isource)%phi)

      call reduce(sources(isource)%theta,sources(isource)%phi)

      ! converts geocentric to geographic colatitude
      call geocentric_2_geographic_dble(sources(isource)%theta,sources(isource)%colatitude)

      ! brings longitude between -PI and PI
      if (sources(isource)%phi > PI)sources(isource)%phi = sources(isource)%phi - TWO_PI

      write(IMAIN,*)
      write(IMAIN,*) '  original (requested) position of the source:'
      write(IMAIN,*)
      write(IMAIN,*) '        latitude: ',sources(isource)%latitude
      write(IMAIN,*) '       longitude: ',sources(isource)%longitude
      write(IMAIN,*) '           depth: ',sources(isource)%depth,' km'
      write(IMAIN,*)

      ! compute real position of the source
      write(IMAIN,*) '  position of the source that will be used:'
      write(IMAIN,*)
      write(IMAIN,*) '        latitude: ',(PI_OVER_TWO-sources(isource)%colatitude)*RADIANS_TO_DEGREES
      write(IMAIN,*) '       longitude: ', sources(isource)%phi*RADIANS_TO_DEGREES
      write(IMAIN,*) '           depth: ',(sources(isource)%r0-sources(isource)%r_found)*R_PLANET/1000.0d0,' km'
      write(IMAIN,*)

      ! display error in location estimate
      write(IMAIN,*) '  Error in location of the source: ',sngl(sources(isource)%final_distance),' km'

      ! add warning if estimate is poor
      ! (usually means source outside the mesh given by the user)
      if (sources(isource)%final_distance > 5.d0) then
        write(IMAIN,*)
        write(IMAIN,*) '*****************************************************'
        write(IMAIN,*) '*****************************************************'
        write(IMAIN,*) '***** WARNING: source location estimate is poor *****'
        write(IMAIN,*) '*****************************************************'
        write(IMAIN,*) '*****************************************************'
      endif
      call flush(IMAIN)


    enddo ! end of loop over all source subsets

    ! display maximum error in location estimate

    ! sets total magnitude (for finite sources)
    M0 = total_M0
    Mw = total_Mw
    force_N = total_force_N

    if (NSOURCES > 1) then
      write(IMAIN,*)
      write(IMAIN,*) '********************'
      write(IMAIN,*) 'finite source combined over all ',NSOURCES,' sources applied:'
      if (sources(isource)%force) then
        ! total force in Newton
        write(IMAIN,*) '  total force = ', sngl(force_N),'(Newton)' ! dimensionalized
      else
        ! moment-tensor
        write(IMAIN,*) '     total scalar moment M0 = ', M0,' dyne-cm'
        write(IMAIN,*) '  total moment magnitude Mw = ', Mw
      endif
      write(IMAIN,*) '********************'
    endif

    ! compute maximal distance for all the sources
    source_final_distance_max = ZERO
    do isource=1,NSOURCES
      if (sources(isource)%final_distance > source_final_distance_max) then
        source_final_distance_max = sources(isource)%final_distance
      endif
    enddo

    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of the sources: ',sngl(source_final_distance_max),' km'
    write(IMAIN,*)



    end subroutine locate_sources


  !
  !-------------------------------------------------------------------------------------------------
  !

  ! subroutine read_source_locations(sources)

  !   ! forces
  !   use ctypes, only: t_source
  !   implicit none

  !   type(t_source), dimension(:), intent(inout) :: sources

  !   ! local parameters
  !   ! event time
  !   integer :: yr,jda,mo,da,ho,mi
  !   double precision :: sec

  !   ! initializes
  !   srclat(:) = 0.d0
  !   srclon(:) = 0.d0
  !   srcdepth(:) = 0.d0
  !   moment_tensor(:,:) = 0.d0

  !   tshift_src(:) = 0.d0
  !   hdur(:) = 0.d0
  !   min_tshift_src_original = 0.d0

  !   ! reads in source descriptions
  !   if (USE_FORCE_POINT_SOURCE) then
  !     ! point forces
  !     if (myrank == 0) then
  !       ! only main process reads in FORCESOLUTION file
  !       call get_force(tshift_src,hdur, &
  !                      srclat,srclon,srcdepth,DT,NSOURCES, &
  !                      min_tshift_src_original,force_stf,factor_force_source, &
  !                      comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP)
  !     endif
  !     ! broadcasts specific point force infos
  !     call bcast_all_i(force_stf,NSOURCES)
  !     call bcast_all_dp(factor_force_source,NSOURCES)
  !     call bcast_all_dp(comp_dir_vect_source_E,NSOURCES)
  !     call bcast_all_dp(comp_dir_vect_source_N,NSOURCES)
  !     call bcast_all_dp(comp_dir_vect_source_Z_UP,NSOURCES)
  !   else
  !     ! CMT moment tensors
  !     if (myrank == 0) then
  !       ! only main process reads in CMTSOLUTION file
  !       call get_cmt(yr,jda,mo,da,ho,mi,sec, &
  !                    tshift_src,hdur, &
  !                    srclat,srclon,srcdepth,moment_tensor, &
  !                    DT,NSOURCES,min_tshift_src_original)
  !     endif
  !     ! broadcast ispecific moment tensor infos
  !     call bcast_all_dp(moment_tensor,6*NSOURCES)
  !   endif

  !   ! broadcast the information read on the main node to all the nodes
  !   call bcast_all_dp(tshift_src,NSOURCES)
  !   call bcast_all_dp(hdur,NSOURCES)
  !   call bcast_all_dp(srclat,NSOURCES)
  !   call bcast_all_dp(srclon,NSOURCES)
  !   call bcast_all_dp(srcdepth,NSOURCES)
  !   call bcast_all_singledp(min_tshift_src_original)

  ! end subroutine read_source_locations

  !
  !-------------------------------------------------------------------------------------------------
  !

end module source_location