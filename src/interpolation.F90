module interpolation

  private
  public :: interpolateMT


contains

  subroutine interpolateMT(&
    displacement, NSTAT, NT, NGLLX, NGLLY, NGLLZ, xigll, yigll, zigll, &
    Mxx, Myy, Mzz, Mxy, Mxz, Myz, xi,eta,gamma, &
    xix, xiy, xiz, etax, etay, etaz, gammax, gammay, gammaz, &
    seismograms)

    use lagrange_poly, only: lagrange_any

    implicit none

    ! IN
    double precision, dimension(:,:,:,:,:,:,:), intent(in) :: displacement
    integer :: NSTAT, NGLLX, NGLLY, NGLLZ
    integer(kind=8) :: NT
    double precision :: Mxx, Myy, Mzz, Mxy, Mxz, Myz
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
    double precision, dimension(NGLLX) :: hxi, hpxi
    double precision, dimension(NGLLY) :: heta, hpeta
    double precision, dimension(NGLLZ) :: hgamma, hpgamma
    double precision :: h_x, h_y, h_z, h_xi, h_eta, h_gamma
    double precision, dimension(6) :: moment_tensor
    double precision, dimension(:,:,:,:), allocatable :: epsilon
    integer :: istat, icomp, im, i, j, k
    integer :: NCOMP = 3

    allocate(epsilon(NSTAT, NCOMP, 6, NT))

    ! Get interpolation values
    call lagrange_any(xi,NGLLX,xigll,hxi,hpxi)
    call lagrange_any(eta,NGLLY,yigll,heta,hpeta)
    call lagrange_any(gamma,NGLLZ,zigll,hgamma,hpgamma)


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

          do istat=1,NSTAT
            do icomp=1,NCOMP
              epsilon(istat, icomp, 1, :) = epsilon(istat, icomp, 1, :) + (displacement(istat,icomp,1,i,j,k,:) * h_x)
              epsilon(istat, icomp, 2, :) = epsilon(istat, icomp, 2, :) + (displacement(istat,icomp,2,i,j,k,:) * h_y)
              epsilon(istat, icomp, 3, :) = epsilon(istat, icomp, 3, :) + (displacement(istat,icomp,3,i,j,k,:) * h_z)
              epsilon(istat, icomp, 4, :) = epsilon(istat, icomp, 4, :) + 0.5 * (displacement(istat,icomp,2,i,j,k,:) * h_x &
                                                    + displacement(istat,icomp,1,i,j,k,:) * h_y)
              epsilon(istat, icomp, 5, :) = epsilon(istat, icomp, 5, :) + 0.5 * (displacement(istat,icomp,3,i,j,k,:) * h_x &
                                                    + displacement(istat,icomp,1,i,j,k,:) * h_z)
              epsilon(istat, icomp, 6, :) = epsilon(istat, icomp, 6, :) + 0.5 * (displacement(istat,icomp,3,i,j,k,:) * h_y &
                                                    + displacement(istat,icomp,2,i,j,k,:) * h_z)
            enddo
          enddo
        enddo
      enddo
    enddo

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


    deallocate(epsilon)

  end subroutine interpolateMT

end module interpolation