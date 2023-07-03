module interpolation

  private
  public :: interpolateMT, spline1d, interp1d, spline1d_seismograms


contains

  subroutine spline1d_seismograms(t, seismograms, tq, seisq)
    implicit none

    ! Arguments
    double precision, dimension(:), intent(in) :: t
    double precision, dimension(:), intent(in) :: tq
    double precision, dimension(:,:,:), intent(in) :: seismograms

    ! OUT
    double precision, dimension(:,:,:), intent(out) :: seisq

    ! Local
    integer :: istat, icomp
    integer :: Nstat, Ncomp

    ! Get number of stations
    Nstat = size(seismograms, dim=1)
    Ncomp = size(seismograms, dim=2)

    ! Interpolate using splines.
    do istat=1,Nstat
      do icomp=1,Ncomp
        call spline1d(t,seismograms(istat, icomp, :),tq, seisq(istat, icomp, :))
      enddo
    enddo


  end subroutine spline1d_seismograms

  ! ============================================================================
  ! ============================================================================
  ! ============================================================================

  subroutine spline1d(x,y,xq,yq)
    use splines, only: spline
    implicit none

    ! Arguments
    double precision, dimension(:), intent(in) :: x, y
    double precision, dimension(:), intent(in) :: xq
    double precision, dimension(:), intent(out) :: yq

    ! Local values
    integer :: ix, Nx
    type(spline) :: sp

    ! Get size of the array
    Nx = size(xq)

    ! Create spline coefficients
    sp = spline(x,y)

    ! Check whether min/max of xq is within x
    if (.false. .eqv. (((minval(x) .le. minval(xq)) .and. (minval(xq) .le. maxval(x))) &
      .and. ((minval(x) .le. maxval(xq)) .and. (maxval(xq) .le. maxval(x))))) then
      stop "Error in spline interpolation. Interpolation range not in data range."
    endif

    ! Interpolate Values.
    do ix=1,Nx
      yq(ix) = sp%value(xq(ix))
    enddo

  end subroutine spline1d

  ! ============================================================================
  ! ============================================================================
  ! ============================================================================

  subroutine interp1d( xData, yData, xVal, yVal )
  !
  ! Linear 1D interpolation.
  !
  ! Inputs: xData = a vector of the x-values of the data to be interpolated
  !         yData = a vector of the y-values of the data to be interpolated
  !         xVal  = a vector of the x-values where interpolation should be performed
  ! Output: yVal  = a vector of the resulting interpolated values
  !
  ! Disclaimers: Only works for monotonically-increasing interpolation values
  !
  !

    implicit none
    ! Args
    double precision, intent(in)  :: xData(:)                 ! Input x-data
    double precision, intent(in)  :: yData(:)                 ! Input y-data
    double precision, intent(in)  :: xVal(:)                  ! Vector with uniformly spaced query values
    double precision, intent(out) :: yVal(:)                  ! Output queried values

    ! Locals
    integer :: inputIndex, dataIndex
    double precision :: minXdata, maxXdata,  minXVal, maxXVal, xRange, weight

    ! Check array sizes
    if (.not. (size(xData) == size(yData))) stop "x and y data have to have the same size."
    if (.not. (size(xVal) == size(yVal)))   stop "x and y query points have to have the same size."

    ! Check if monotonically increasing
    if (minval(xData(2:size(xData)) - xData(1:size(xData)-1)) .le. 0.d0) then
      stop "Data not monotonically increasing."
    endif

    ! Check if monotonically increasing
    if (minval(xVal(2:size(xVal)) - xVal(1:size(xVal)-1)) .le. 0.d0) then
      stop "x query values not monotonically increasing."
    endif

    ! Check range
    minXData = xData(1)
    maxXData = xData(size(xData))
    minXVal = xVal(1)
    maxXVal = xVal(size(xVal))

    ! Check whether min/max of xq is within x
    if (.false. .eqv. (((minXData .le. minXVal) .and. (minXVal .le. maxXData)) &
      .and. ((minXData .le. maxXVal) .and. (maxXVal .le. maxXData)))) then
      stop "Error in spline interpolation. Interpolation range not in data range."
    endif


    do inputIndex = 1, size(xVal)
      ! possible checks for out of range xVal could go here

      ! this will work if x is uniformly spaced, otherwise increment
      ! dataIndex until xData(dataIndex+1)>xVal(inputIndex)
      dataIndex = 1
      do while (xData(dataIndex+1) < xVal(inputIndex))
        dataIndex = dataIndex + 1
      enddo

      weight = (xVal(inputIndex) - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex))
      yVal(inputIndex) = (1.0-weight)*yData(dataIndex) + weight*yData(dataIndex+1)
    enddo
  end subroutine interp1d

  ! ============================================================================
  ! ============================================================================
  ! ============================================================================

  subroutine interpolateMT(&
    displacement, NSTAT, NT, NGLLX, NGLLY, NGLLZ, xigll, yigll, zigll, &
    Mxx, Myy, Mzz, Mxy, Mxz, Myz, xi,eta,gamma, &
    xix, xiy, xiz, etax, etay, etaz, gammax, gammay, gammaz, &
    seismograms)

    use lagrange_poly, only: lagrange_any
    use constants, only: DEBUG
    implicit none

    ! IN
    double precision, dimension(:,:,:,:,:,:,:), intent(in) :: displacement
    integer :: NSTAT, NGLLX, NGLLY, NGLLZ
    integer(kind=8) :: NT
    real :: start, finish
    real :: startsub, finishsub
    double precision, intent(in) :: Mxx, Myy, Mzz, Mxy, Mxz, Myz
    double precision,intent(in) :: xi,eta,gamma
    double precision, intent(in) :: xix, xiy, xiz
    double precision, intent(in) :: etax, etay, etaz
    double precision, intent(in) :: gammax, gammay, gammaz
    double precision, dimension(NGLLX) :: xigll
    double precision, dimension(NGLLY) :: yigll
    double precision, dimension(NGLLZ) :: zigll

    ! OUT
    double precision, dimension(:,:,:), intent(out) :: seismograms

    ! Local
    integer,parameter :: NCOMP = 3
    double precision, dimension(NGLLX) :: hxi, hpxi
    double precision, dimension(NGLLY) :: heta, hpeta
    double precision, dimension(NGLLZ) :: hgamma, hpgamma
    double precision :: h_x, h_y, h_z, h_xi, h_eta, h_gamma
    double precision, dimension(6) :: moment_tensor
    double precision, dimension(NSTAT, NCOMP, 6, NT) :: epsilon
    integer :: istat, icomp, im, i, j, k

    if (DEBUG) call cpu_time(start)

    ! Initialize(!) we are summing over this one!
    epsilon(:,:,:,:) = 0.d0


    ! Get interpolation values
    call lagrange_any(xi,NGLLX,xigll,hxi,hpxi)
    call lagrange_any(eta,NGLLY,yigll,heta,hpeta)
    call lagrange_any(gamma,NGLLZ,zigll,hgamma,hpgamma)


    if (DEBUG) call cpu_time(startsub)

    do k=1, NGLLZ
      do j=1, NGLLY
        do i=1, NGLLX

          ! Compute derivative with respect to element coordinates
          h_xi = hpxi(i)* heta(j) * hgamma(k)
          h_eta = hxi(i) * hpeta(j) * hgamma(k)
          h_gamma = hxi(i) * heta(j) * hpgamma(k)

          ! Transform derivative using the Jacobian values
          h_x = h_xi * xix + h_eta * etax + h_gamma * gammax
          h_y = h_xi * xiy + h_eta * etay + h_gamma * gammay
          h_z = h_xi * xiz + h_eta * etaz + h_gamma * gammaz

          ! do istat=1,NSTAT
          !   do icomp=1,NCOMP
          epsilon(:, :, 1, :) = epsilon(:, :, 1, :) + (displacement(:,:,1,i,j,k,:) * h_x)
          epsilon(:, :, 2, :) = epsilon(:, :, 2, :) + (displacement(:,:,2,i,j,k,:) * h_y)
          epsilon(:, :, 3, :) = epsilon(:, :, 3, :) + (displacement(:,:,3,i,j,k,:) * h_z)
          epsilon(:, :, 4, :) = epsilon(:, :, 4, :) + 0.5 * (displacement(:,:,2,i,j,k,:) * h_x &
                                                + displacement(:,:,1,i,j,k,:) * h_y)
          epsilon(:, :, 5, :) = epsilon(:, :, 5, :) + 0.5 * (displacement(:,:,3,i,j,k,:) * h_x &
                                                + displacement(:,:,1,i,j,k,:) * h_z)
          epsilon(:, :, 6, :) = epsilon(:, :, 6, :) + 0.5 * (displacement(:,:,3,i,j,k,:) * h_y &
                                                + displacement(:,:,2,i,j,k,:) * h_z)
          !   enddo
          ! enddo
        enddo
      enddo
    enddo

    if (DEBUG) call cpu_time(finishsub)
    if (DEBUG) print '("    ____mainloop took ",f6.3," seconds.")', finishsub-startsub


    ! Initialize dot product
    moment_tensor(1) = Mxx
    moment_tensor(2) = Myy
    moment_tensor(3) = Mzz
    moment_tensor(4) = 2*Mxy
    moment_tensor(5) = 2*Mxz
    moment_tensor(6) = 2*Myz

    ! Initialize seismograms
    seismograms(:,:,:) = 0.d0
    do im=1,6
      seismograms(:,:,:) = seismograms(:,:,:) + epsilon(:,:,im,:) * moment_tensor(im)
    enddo

    if (DEBUG) call cpu_time(finish)
    if (DEBUG) print '("    interpolateMT took ",f6.3," seconds.")', finish-start

  end subroutine interpolateMT

end module interpolation