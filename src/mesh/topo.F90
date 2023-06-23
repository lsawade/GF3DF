
module topo

  implicit none
  private
  public :: get_topo_bathy


contains
  subroutine get_topo_bathy(xlat,xlon,value, NX_BATHY, NY_BATHY, ibathy_topo, RESOLUTION_TOPO_FILE)

  !
  !---- get elevation or ocean depth in meters at a given latitude and longitude
  !
    use constants, only: ZERO

    implicit none

    ! location latitude/longitude (in degree)
    double precision,intent(in):: xlat,xlon

    ! returns elevation (in meters)
    double precision,intent(out):: value

    ! use integer array to store values
    integer, intent(in) :: NX_BATHY, NY_BATHY
    integer, dimension(:,:),intent(in) :: ibathy_topo

    real(kind=8), intent(in) :: RESOLUTION_TOPO_FILE

    ! local parameters
    integer :: iadd1,iel1

    double precision :: samples_per_degree_topo
    double precision :: xlo
    double precision :: lon_corner,lat_corner,ratio_lon,ratio_lat

    ! initializes elevation
    value = ZERO

    ! longitude within range [0,360] degrees
    xlo = xlon
    if (xlo < 0.d0) xlo = xlo + 360.d0
    if (xlo > 360.d0) xlo = xlo - 360.d0

    ! compute number of samples per degree
    samples_per_degree_topo = dble(RESOLUTION_TOPO_FILE) / 60.d0

    ! compute offset in data file and avoid edge effects
    iadd1 = 1 + int((90.d0-xlat)/samples_per_degree_topo)
    if (iadd1 < 1) iadd1 = 1
    if (iadd1 > NY_BATHY) iadd1 = NY_BATHY

    iel1 = int(xlo/samples_per_degree_topo)
    if (iel1 <= 0 .or. iel1 > NX_BATHY) iel1 = NX_BATHY

  ! Use bilinear interpolation rather nearest point interpolation

    ! convert integer value to double precision
    !  value = dble(ibathy_topo(iel1,iadd1))

    lon_corner = iel1 * samples_per_degree_topo
    lat_corner = 90.d0 - iadd1 * samples_per_degree_topo

    ratio_lon = (xlo-lon_corner)/samples_per_degree_topo
    ratio_lat = (xlat-lat_corner)/samples_per_degree_topo

    if (ratio_lon < 0.d0) ratio_lon = 0.d0
    if (ratio_lon > 1.d0) ratio_lon = 1.d0
    if (ratio_lat < 0.d0) ratio_lat = 0.d0
    if (ratio_lat > 1.d0) ratio_lat = 1.d0

    ! convert integer value to double precision
    if (iadd1 <= NY_BATHY-1 .and. iel1 <= NX_BATHY-1) then
      ! interpolates for points within boundaries
      value = dble(ibathy_topo(iel1,iadd1))     * (1.d0-ratio_lon) * (1.d0-ratio_lat) &
            + dble(ibathy_topo(iel1+1,iadd1))   * ratio_lon * (1.d0-ratio_lat) &
            + dble(ibathy_topo(iel1+1,iadd1+1)) * ratio_lon * ratio_lat &
            + dble(ibathy_topo(iel1,iadd1+1))   * (1.d0-ratio_lon) * ratio_lat

    else if (iadd1 <= NY_BATHY-1 .and. iel1 == NX_BATHY) then
      ! interpolates for points on longitude border
      value = dble(ibathy_topo(iel1,iadd1))   * (1.d0-ratio_lon)*(1.d0-ratio_lat) &
            + dble(ibathy_topo(1,iadd1))      * ratio_lon*(1.d0-ratio_lat) &
            + dble(ibathy_topo(1,iadd1+1))    * ratio_lon*ratio_lat &
            + dble(ibathy_topo(iel1,iadd1+1)) * (1.d0-ratio_lon)*ratio_lat

    else
      ! for points on latitude boundaries
      value = dble(ibathy_topo(iel1,iadd1))
    endif

  end subroutine get_topo_bathy

end module topo