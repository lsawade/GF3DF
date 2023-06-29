!
! Standalone program to compute the ellipticity splines. To showcase both eta
! and epsilon
!


  module par

    double precision, parameter :: HOURS_PER_DAY = 24
    double precision , parameter:: SECONDS_PER_HOUR = 3600

    logical :: ONE_CRUST = .true.
    double precision, parameter :: EARTH_RHOAV = 5514.3
    double precision, parameter :: EARTH_R = 6371000.0
    double precision, parameter :: R_PLANET_KM = EARTH_R/1000.d0
    double precision, parameter :: R_PLANET = EARTH_R

    double precision, parameter :: RHOAV = EARTH_RHOAV
    double precision, parameter :: GRAV = 6.67384e-11

    double precision, parameter :: PI = 3.141592653589793d0
    double precision, parameter :: TWO_PI = 2.d0 * PI
    double precision, parameter :: R_UNIT_SPHERE = 1.d0


    ! Number of sample points
    integer, parameter :: NR_DENSITY = 640

    double precision, parameter :: PREM_ROCEAN = 6368000.d0           ! at   3 km depth
    double precision, parameter :: PREM_RSURFACE = 6371000.0
    double precision, parameter :: PREM_RMIDDLE_CRUST = 6356000.0  ! 15 km depth
    double precision, parameter :: PREM_RMOHO = PREM_RSURFACE - 24400.0  ! 6346600.  at 24.4km depth
    double precision, parameter :: PREM_R80 = 6291000.0  ! 80 km
    double precision, parameter :: PREM_R220 = 6151000.0
    double precision, parameter :: PREM_R400 = 5971000.0
    double precision, parameter :: PREM_R600 = 5771000.0
    double precision, parameter :: PREM_R670 = 5701000.0                 ! at 670 km depth
    double precision, parameter :: PREM_R771 = 5600000.0
    double precision, parameter :: PREM_RTOPDDOUBLEPRIME = 3630000.0
    double precision, parameter :: PREM_RCMB = 3480000.0  ! 2891 km depth
    ! note: SPECFEM versions up to 8.0 (Aug, 2021) used an inner core radius of 1221 km
    !       based on Table 1 in Dziewonski & Anderson's PREM paper, the inner core radius is at 1221.5 km
    !       double precision, parameter: : PREM_RICB = 1221000.            ! old versions
    double precision, parameter :: PREM_RICB = 1221500.0

    ! PREM2 additional radii for modifications
    double precision, parameter :: PREM2_RDDOUBLEPRIME_UPPER = 3840000.0  ! upper D'' region at 3840km radius
    double precision, parameter :: PREM2_ROC_LOWER = 1621500.0            ! lower outer core at 1621.5km radius
    double precision, parameter :: PREM2_RIC_UPPER = 1010000.0            ! upper inner core at 1010 km radius

  end module par




program make_ellipticity_prog


  use par, only: NR_DENSITY

  integer :: nspl
  double precision,dimension(NR_DENSITY) :: rspl,ellipicity_spline,ellipicity_spline2

  call make_ellipticity(nspl,rspl,ellipicity_spline,ellipicity_spline2)


contains


  subroutine make_ellipticity(nspl,rspl,ellipicity_spline,ellipicity_spline2)

! creates a spline for the ellipticity profile in PREM
! radius and density are non-dimensional
! or
! for mars, creates a spline for the ellipticity profile in Mars (Sohl $ Spohn)


!   use constants, only: NR_DENSITY,TWO_PI,PI,GRAV,R_UNIT_SPHERE,myrank

!   use shared_parameters, only: PLANET_TYPE,IPLANET_EARTH,IPLANET_MARS,IPLANET_MOON, &
!     R_PLANET,R_PLANET_KM,RHOAV, &
!     HOURS_PER_DAY,SECONDS_PER_HOUR

  use par

  implicit none




  ! reference models
  ! use model_prem_par
  ! use model_sohl_par
  ! use model_vpremoon_par

  integer,intent(inout) :: nspl
  double precision,dimension(NR_DENSITY),intent(inout) :: rspl,ellipicity_spline,ellipicity_spline2

  ! local parameters
  integer :: i
  ! radii
  double precision :: ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670, &
                      R771,RTOPDDOUBLEPRIME,RCMB,RICB,RSURFACE
  double precision :: r_icb,r_cmb,r_topddoubleprime,r_771,r_670,r_600
  double precision :: r_400,r_220,r_80,r_moho,r_middle_crust,r_ocean,r_0
  double precision :: SOHL_RMOHO,SOHL_R80,SOHL_R220,SOHL_R400,SOHL_R600,SOHL_R670, &
                      SOHL_R771,SOHL_RTOPDDOUBLEPRIME,SOHL_RCMB

  double precision,dimension(NR_DENSITY) :: r,rho,epsilonval,eta
  double precision,dimension(NR_DENSITY) :: radau,k
  double precision,dimension(NR_DENSITY) :: s1,s2,s3

  double precision :: z,g_a,bom,exponentval,integral_rho,integral_radau
  double precision :: yp1,ypn

  ! debugging
  logical, parameter :: DEBUG = .true.
  double precision :: ell,radius

! note: ellipicity is implemented by Radau's approximation to Clairaut's equation.
!
!       details can be found in: doc/notes/ellipticity_equations_from_Dahlen_Tromp_1998.pdf
!
! for calculating the ellipicity factor epsilon(r) anywhere on Earth, one needs to integrate over Earth's density profile.
! we use PREM's density model for that (including the ocean layer on top).
!
! as a todo in future: this PREM density profile might be slightly off for different Earth models,
!                      please check the effect.

  ! Earth
  ! PREM radius of the Earth for gravity calculation
  double precision, parameter :: R_EARTH_ELLIPTICITY = 6371000.d0
  ! PREM radius of the ocean floor for gravity calculation
  double precision, parameter :: ROCEAN_ELLIPTICITY = PREM_ROCEAN

  ! Earth
  ! default PREM
  ! radius of the planet for gravity calculation
  RSURFACE = R_EARTH_ELLIPTICITY  ! physical surface (Earth: 6371000, ..)
  ROCEAN = ROCEAN_ELLIPTICITY
  RMIDDLE_CRUST = PREM_RMIDDLE_CRUST
  RMOHO = PREM_RMOHO
  R80  = PREM_R80
  R220 = PREM_R220
  R400 = PREM_R400
  R600 = PREM_R600
  R670 = PREM_R670
  R771 = PREM_R771
  RTOPDDOUBLEPRIME = PREM_RTOPDDOUBLEPRIME
  RCMB = PREM_RCMB
  RICB = PREM_RICB

  ! non-dimensionalize
  r_icb = RICB / RSURFACE
  r_cmb = RCMB / RSURFACE
  r_topddoubleprime = RTOPDDOUBLEPRIME / RSURFACE
  r_771 = R771 / RSURFACE
  r_670 = R670 / RSURFACE
  r_600 = R600 / RSURFACE
  r_400 = R400 / RSURFACE
  r_220 = R220 / RSURFACE
  r_80 = R80 / RSURFACE
  r_moho = RMOHO / RSURFACE
  r_middle_crust = RMIDDLE_CRUST / RSURFACE
  r_ocean = ROCEAN / RSURFACE
  r_0 = 1.d0

  ! sets sampling points in different layers
  ! inner core
  do i = 1,163
    r(i) = r_icb*dble(i-1)/dble(162)
  enddo
  ! outer core
  do i = 164,323
    r(i) = r_icb+(r_cmb-r_icb)*dble(i-164)/dble(159)
  enddo
  ! D''
  do i = 324,336
    r(i) = r_cmb+(r_topddoubleprime-r_cmb)*dble(i-324)/dble(12)
  enddo
  ! D'' to 771
  do i = 337,517
    r(i) = r_topddoubleprime+(r_771-r_topddoubleprime)*dble(i-337)/dble(180)
  enddo
  ! 771 to 670
  do i = 518,530
    r(i) = r_771+(r_670-r_771)*dble(i-518)/dble(12)
  enddo
  ! 670 to 600
  do i = 531,540
    r(i) = r_670+(r_600-r_670)*dble(i-531)/dble(9)
  enddo
  ! 600 to 400
  do i = 541,565
    r(i) = r_600+(r_400-r_600)*dble(i-541)/dble(24)
  enddo
  ! 400 to 220
  do i = 566,590
    r(i) = r_400+(r_220-r_400)*dble(i-566)/dble(24)
  enddo
  ! 220 to 80
  do i = 591,609
    r(i) = r_220+(r_80-r_220)*dble(i-591)/dble(18)
  enddo
  ! 80 to Moho
  do i = 610,619
    r(i) = r_80+(r_moho-r_80)*dble(i-610)/dble(9)
  enddo
  ! Moho to middle crust
  do i = 620,626
    r(i) = r_moho+(r_middle_crust-r_moho)*dble(i-620)/dble(6)
  enddo
  ! middle crust to ocean
  do i = 627,633
    r(i) = r_middle_crust+(r_ocean-r_middle_crust)*dble(i-627)/dble(6)
  enddo
  ! ocean
  do i = 634,NR_DENSITY   ! NR_DENSITY = 640
    r(i) = r_ocean+(r_0-r_ocean)*dble(i-634)/dble(6)
  enddo

  ! Earth
  ! use PREM to get the density profile for ellipticity (fine for other 1D reference models)
  do i = 1,NR_DENSITY
    call prem_density(r(i),rho(i))
    radau(i) = rho(i)*r(i)*r(i)
  enddo

  eta(1) = 0.0d0
  k(1) = 0.0d0

  do i = 2,NR_DENSITY
    call intgrl(integral_rho,r,1,i,rho,s1,s2,s3)

! Radau approximation of Clairaut's equation for first-order terms of ellipticity, see e.g. Jeffreys H.,
! The figures of rotating planets, Mon. Not. R. astr. Soc., vol. 113, p. 97-105 (1953).
! The Radau approximation is mentioned on page 97.
! For more details see Section 14.1.2 in Dahlen and Tromp (1998)
! (see also in file ellipticity_equations_from_Dahlen_Tromp_1998.pdf in the "doc" directory of the code).
    call intgrl(integral_radau,r,1,i,radau,s1,s2,s3)

    z = (2.0d0/3.0d0) * integral_radau / (integral_rho*r(i)*r(i))

    ! this comes from equation (14.19) in Dahlen and Tromp (1998)
    eta(i) = (25.0d0/4.0d0)*((1.0d0-(3.0d0/2.0d0)*z)**2.0d0)-1.0d0
    k(i) = eta(i)/(r(i)**3.0d0)
  enddo

  ! day rotation
  bom = TWO_PI/(HOURS_PER_DAY*SECONDS_PER_HOUR)

  ! non-dimensionalized value
  bom = bom/sqrt(PI*GRAV*RHOAV)

  ! Integral rho must be including PI??
  g_a = 4.0d0 * integral_rho

  ! this is the equation right above (14.21) in Dahlen and Tromp (1998)
  epsilonval(NR_DENSITY) = (5.0d0/2.d0)*(bom**2.0d0)*R_UNIT_SPHERE / (g_a * (eta(NR_DENSITY)+2.0d0))

  do i = 1,NR_DENSITY-1
    call intgrl(exponentval,r,i,NR_DENSITY,k,s1,s2,s3)
    epsilonval(i) = epsilonval(NR_DENSITY)*exp(-exponentval)
  enddo

  ! initializes spline coefficients
  rspl(:) = 0.d0
  ellipicity_spline(:) = 0.d0
  ellipicity_spline2(:) = 0.d0

  ! get ready to spline epsilonval
  nspl = 1
  rspl(1) = r(1)
  ellipicity_spline(1) = epsilonval(1)
  do i = 2,NR_DENSITY
    if (r(i) /= r(i-1)) then
      nspl = nspl+1
      rspl(nspl) = r(i)
      ellipicity_spline(nspl) = epsilonval(i)
    endif
  enddo

  ! spline epsilonval
  yp1 = 0.0d0
  ypn = (5.0d0/2.0d0)*(bom**2)/g_a-2.0d0*epsilonval(NR_DENSITY)

  call spline_construction(rspl,ellipicity_spline,nspl,yp1,ypn,ellipicity_spline2)

  ! debug
  if (DEBUG) then

    ! print *,'debug: make ellipticity'
    ! print *,'debug: number of splines = ',nspl
    print *,'#r(km) #ellipticity_factor #eta #radau #k #r(non-dim) #i'
    do i = 1,NR_DENSITY
      radius = r(i)
      ! get ellipticity using spline evaluation
      call spline_evaluation(rspl,ellipicity_spline,ellipicity_spline2,nspl,radius,ell)
      print *,radius*R_PLANET_KM,ell,eta(i), radau(i),k(i),radius,i
    enddo
    print *

  endif

  end subroutine make_ellipticity

  ! =========================================================================
  ! =========================================================================
  ! =========================================================================
  ! =========================================================================
  ! =========================================================================


  subroutine prem_density(x,rho)

  use par

  implicit none

  double precision,intent(in) :: x
  double precision,intent(out) :: rho

  ! local parameters
  double precision :: r

  ! compute real physical radius in meters
  r = x * R_PLANET

  ! calculates density according to radius
  if (r <= PREM_RICB) then
    rho = 13.0885d0 - 8.8381d0*x*x
  else if (r > PREM_RICB .and. r <= PREM_RCMB) then
    rho = 12.5815d0 - 1.2638d0*x - 3.6426d0*x*x - 5.5281d0*x*x*x
  else if (r > PREM_RCMB .and. r <= PREM_RTOPDDOUBLEPRIME) then
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
  else if (r > PREM_RTOPDDOUBLEPRIME .and. r <= PREM_R771) then
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
  else if (r > PREM_R771 .and. r <= PREM_R670) then
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
  else if (r > PREM_R670 .and. r <= PREM_R600) then
    rho = 5.3197d0 - 1.4836d0*x
  else if (r > PREM_R600 .and. r <= PREM_R400) then
    rho = 11.2494d0 - 8.0298d0*x
  else if (r > PREM_R400 .and. r <= PREM_R220) then
    rho = 7.1089d0 - 3.8045d0*x
  else if (r > PREM_R220 .and. r <= PREM_R80) then
    rho = 2.6910d0 + 0.6924d0*x
  else
    if (r > PREM_R80 .and. r <= PREM_RMOHO) then
      rho = 2.6910d0 + 0.6924d0*x
    else if (r > PREM_RMOHO .and. r <= PREM_RMIDDLE_CRUST) then
      if (ONE_CRUST) then
        rho = 2.6d0  ! takes upper crust value
      else
        rho = 2.9d0
      endif
    else if (r > PREM_RMIDDLE_CRUST .and. r <= PREM_ROCEAN) then
      rho = 2.6d0
    else if (r > PREM_ROCEAN) then
      rho = 2.6d0 ! extends upper crust
    endif
  endif

  ! non-dimensionalizes
  rho = rho * 1000.0d0 / RHOAV

  end subroutine prem_density


  ! =========================================================================
  ! =========================================================================
  ! =========================================================================
  ! =========================================================================
  ! =========================================================================

  subroutine intgrl(sumval,r,nir,ner,f,s1,s2,s3)

! Computes the integral of f[i]*r[i]*r[i] from i=nir to i=ner for
! radii values as in model PREM_an640

  use par, only: NR_DENSITY
  implicit none

! Argument variables
  integer,intent(in) :: ner,nir
  double precision,dimension(NR_DENSITY),intent(in) :: f,r
  double precision,dimension(NR_DENSITY),intent(inout) :: s1,s2,s3
  double precision,intent(out) :: sumval

! Local variables
  double precision, parameter :: third = 1.0d0/3.0d0
  double precision, parameter :: fifth = 1.0d0/5.0d0
  double precision, parameter :: sixth = 1.0d0/6.0d0

  double precision :: rji,yprime(640)
  double precision :: s1l,s2l,s3l

  integer :: i,j,n,kdis(28)
  integer :: ndis,nir1

  data kdis/163,323,336,517,530,540,565,590,609,619,626,633,16*0/

  ndis = 12
  n = NR_DENSITY

  call deriv(f,yprime,n,r,ndis,kdis,s1,s2,s3)

  nir1 = nir + 1
  sumval = 0.0d0

  do i = nir1,ner
    j = i-1
    rji = r(i) - r(j)
    s1l = s1(j)
    s2l = s2(j)
    s3l = s3(j)
    sumval = sumval &
              + r(j) * r(j) * rji * ( f(j) + rji*(0.5d0*s1l + rji*(third*s2l + rji*0.25d0*s3l)) ) &
              + 2.0d0 * r(j) * rji * rji * ( 0.5d0*f(j) + rji*(third*s1l + rji*(0.25d0*s2l + rji*fifth*s3l)) ) &
              + rji * rji * rji * ( third*f(j) + rji*(0.25d0*s1l + rji*(fifth*s2l + rji*sixth*s3l)) )
  enddo

  end subroutine intgrl

! -------------------------------

  subroutine deriv(y,yprime,n,r,ndis,kdis,s1,s2,s3)

  implicit none

! Argument variables
  integer,intent(in) :: kdis(28),n,ndis
  double precision,intent(in) :: r(n)
  double precision,intent(inout) :: s1(n),s2(n),s3(n)
  double precision,intent(in) :: y(n)
  double precision,intent(inout) :: yprime(n)

! Local variables
  integer :: i,j,j1,j2
  integer :: k,nd,ndp
  double precision :: a0,b0,b1
  double precision :: f(3,1000),h,h2,h2a
  double precision :: h2b,h3a,ha,s13
  double precision :: s21,s32,yy(3)

  yy(1) = 0.d0
  yy(2) = 0.d0
  yy(3) = 0.d0

  ndp = ndis+1
  do 3 nd = 1,ndp
    if (nd == 1) then
      j1 = 1
      j2 = kdis(1)-2
    else if (nd == ndp) then
      j1 = kdis(ndis)+1
      j2 = n-2
    else
      j1 = kdis(nd-1)+1
      j2 = kdis(nd)-2
    endif

    if ((j2+1-j1) > 0) goto 11

    j2 = j2+2
    if (abs(r(j2)-r(j1)) > 0.d0) then
      yy(1) = (y(j2)-y(j1)) / (r(j2)-r(j1))
    else
      yy(1) = 0.d0
    endif
    s1(j1) = yy(1)
    s1(j2) = yy(1)
    s2(j1) = yy(2)
    s2(j2) = yy(2)
    s3(j1) = yy(3)
    s3(j2) = yy(3)
  goto 3

  11 a0 = 0.0d0

    if (j1 == 1) then
      b0 = 0.0d0
    else
      h = r(j1+1)-r(j1)
      h2 = r(j1+2)-r(j1)
      yy(1) = h*h2*(h2-h)
      h = h*h
      h2 = h2*h2
      if (abs(yy(1)) > 0.d0) then
        b0 = (y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/yy(1)
      else
        b0 = 0.d0
      endif
    endif
    b1 = b0

    if (j2 > 1000) stop 'Error in subroutine deriv for j2'

    do i = j1,j2
      h = r(i+1)-r(i)
      yy(1) = y(i+1)-y(i)
      h2 = h*h
      ha = h-a0
      h2a = h-2.0d0*a0
      h3a = 2.0d0*h-3.0d0*a0
      h2b = h2*b0
      s1(i) = h2/ha
      s2(i) = -ha/(h2a*h2)
      s3(i) = -h*h2a/h3a
      f(1,i) = (yy(1)-h*b0)/(h*ha)
      f(2,i) = (h2b-yy(1)*(2.0d0*h-a0))/(h*h2*h2a)
      f(3,i) = -(h2b-3.0d0*yy(1)*ha)/(h*h3a)
      a0 = s3(i)
      b0 = f(3,i)
    enddo

    i = j2+1
    h = r(i+1)-r(i)
    yy(1) = y(i+1)-y(i)
    h2 = h*h
    ha = h-a0
    h2a = h*ha
    h2b = h2*b0-yy(1)*(2.d0*h-a0)
    s1(i) = h2/ha
    f(1,i) = (yy(1)-h*b0)/h2a
    ha = r(j2)-r(i+1)
    yy(1) = -h*ha*(ha+h)
    ha = ha*ha
    yy(1) = (y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/yy(1)
    s3(i) = (yy(1)*h2a+h2b)/(h*h2*(h-2.0d0*a0))
    s13 = s1(i)*s3(i)
    s2(i) = f(1,i)-s13

    do j = j1,j2
      k = i-1
      s32 = s3(k)*s2(i)
      s1(i) = f(3,k)-s32
      s21 = s2(k)*s1(i)
      s3(k) = f(2,k)-s21
      s13 = s1(k)*s3(k)
      s2(k) = f(1,k)-s13
      i = k
    enddo

    s1(i) = b1
    j2 = j2+2
    s1(j2) = yy(1)
    s2(j2) = yy(2)
    s3(j2) = yy(3)
 3 continue

  do i = 1,n
    yprime(i) = s1(i)
  enddo

  end subroutine deriv


  ! =========================================================================
  ! =========================================================================
  ! =========================================================================
  ! =========================================================================
  ! =========================================================================

  subroutine spline_construction(xpoint,ypoint,npoint,tangent_first_point,tangent_last_point,spline_coefficients)

  implicit none

! number of input points and coordinates of the input points
  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(in) :: xpoint,ypoint

! tangent to the spline imposed at the first and last points
  double precision, intent(in) :: tangent_first_point,tangent_last_point

! spline coefficients output by the routine
  double precision, dimension(npoint), intent(out) :: spline_coefficients

  ! local parameters
  integer :: i
  double precision, dimension(:), allocatable :: temporary_array

  allocate(temporary_array(npoint))

  spline_coefficients(1) = - 1.d0 / 2.d0

  temporary_array(1) = (3.d0/(xpoint(2)-xpoint(1)))*((ypoint(2)-ypoint(1))/(xpoint(2)-xpoint(1))-tangent_first_point)

  do i = 2,npoint-1

    spline_coefficients(i) = ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))-1.d0) &
       / ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*spline_coefficients(i-1)+2.d0)

    temporary_array(i) = (6.d0*((ypoint(i+1)-ypoint(i))/(xpoint(i+1)-xpoint(i)) &
       - (ypoint(i)-ypoint(i-1))/(xpoint(i)-xpoint(i-1)))/(xpoint(i+1)-xpoint(i-1)) &
       - (xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*temporary_array(i-1)) &
       / ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*spline_coefficients(i-1)+2.d0)

  enddo

  spline_coefficients(npoint) = ((3.d0/(xpoint(npoint)-xpoint(npoint-1))) &
      * (tangent_last_point-(ypoint(npoint)-ypoint(npoint-1))/(xpoint(npoint)-xpoint(npoint-1))) &
      - 1.d0/2.d0*temporary_array(npoint-1))/(1.d0/2.d0*spline_coefficients(npoint-1)+1.d0)

  do i = npoint-1,1,-1
    spline_coefficients(i) = spline_coefficients(i)*spline_coefficients(i+1) + temporary_array(i)
  enddo

  deallocate(temporary_array)

  end subroutine spline_construction

! --------------

! evaluate a spline

  subroutine spline_evaluation(xpoint,ypoint,spline_coefficients,npoint,x_evaluate_spline,y_spline_obtained)

  use par, only: R_PLANET_KM

  implicit none

! number of input points and coordinates of the input points
  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(in) :: xpoint,ypoint

! spline coefficients to use
  double precision, dimension(npoint), intent(in) :: spline_coefficients

! abscissa at which we need to evaluate the value of the spline
  double precision, intent(in):: x_evaluate_spline

! ordinate evaluated by the routine for the spline at this abscissa
  double precision, intent(out):: y_spline_obtained

  ! local parameters
  integer :: index_loop,index_lower,index_higher
  double precision :: coef1,coef2

! initialize to the whole interval
  index_lower = 1
  index_higher = npoint

! determine the right interval to use, by dichotomy
  do while (index_higher - index_lower > 1)
! compute the middle of the interval
    index_loop = (index_higher + index_lower) / 2
    if (index_loop < 1) index_loop = 1

    if (xpoint(index_loop) > x_evaluate_spline) then
      index_higher = index_loop
    else
      index_lower = index_loop
    endif
  enddo

! test that the interval obtained does not have a size of zero
! (this could happen for instance in the case of duplicates in the input list of points)
  if (xpoint(index_higher) == xpoint(index_lower)) then
    print *,'Error: invalid spline evalution index_higher/index_lower = ',index_higher,index_lower,'range = ',1,npoint
    print *,'       point x = ',x_evaluate_spline,' x radius = ',x_evaluate_spline * R_PLANET_KM,'(km)'
    stop 'incorrect interval found in spline evaluation'
  endif

  coef1 = (xpoint(index_higher) - x_evaluate_spline) / (xpoint(index_higher) - xpoint(index_lower))
  coef2 = (x_evaluate_spline - xpoint(index_lower)) / (xpoint(index_higher) - xpoint(index_lower))

  y_spline_obtained = coef1*ypoint(index_lower) + coef2*ypoint(index_higher) + &
        ((coef1**3 - coef1)*spline_coefficients(index_lower) + &
         (coef2**3 - coef2)*spline_coefficients(index_higher))*((xpoint(index_higher) - xpoint(index_lower))**2)/6.d0

  end subroutine spline_evaluation





end program make_ellipticity_prog