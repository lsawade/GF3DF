submodule (gf3d:seismograms) source_interpolation_fast

  implicit none

contains

  module subroutine interpolate_source_slow(GF, source, seismograms)

    use gf, only: t_GF
    use sources, only: t_source
    use interpolation, only: interpolateMT
    use constants, only: DEBUG

    ! In
    type(t_GF), intent(in) :: GF
    type(t_source), intent(in) :: source
    ! Out
    double precision, dimension(:,:,:) :: seismograms

    ! Local
    double precision, dimension(size(GF%displacement,1), 3, 3, GF%ngllx,GF%nglly,GF%ngllz,GF%nsteps) :: displacement
    integer :: i,j,k,iglob
    real :: start, finish

    if (DEBUG) call cpu_time(start)

    ! Get displacement array
    do k=1,GF%ngllz
      do j=1,GF%nglly
        do i=1,GF%ngllx
          iglob = GF%ibool(i,j,k, source%ispec)
          displacement(:,:,:,i,j,k, :) = GF%displacement(:,:,:,iglob,:)
        enddo
      enddo
    enddo

    if (source%force) then
      call throwerror(-1, "Force source interpolation not implemented.")
    else
      call interpolateMT(&
          displacement, size(GF%displacement,1), &
          GF%nsteps, GF%ngllx, GF%nglly, GF%ngllz, GF%xigll, GF%yigll, GF%zigll, &
          source%Mxx, source%Myy, source%Mzz, &
          source%Mxy, source%Mxz, source%Myz, &
          source%xi, source%eta,source%gamma, &
          source%xix, source%xiy, source%xiz, &
          source%etax, source%etay, source%etaz, &
          source%gammax, source%gammay, source%gammaz, &
          seismograms)
    endif

    if (DEBUG) call cpu_time(finish)
    if (DEBUG) print '("interpolate_source took ",f6.3," seconds.")', finish-start

  end subroutine interpolate_source_slow


  module subroutine interpolate_source(GF, source, seismograms)

    use gf, only: t_GF
    use sources, only: t_source
    use interpolation, only: interpolateMT
    use constants, only: DEBUG
    use lagrange_poly, only: lagrange_any

    implicit none

    ! In
    type(t_GF), intent(in) :: GF
    type(t_source), intent(in) :: source
    ! Out
    double precision, dimension(:,:,:) :: seismograms

    ! Local
    integer :: i,j,k,iglob
    integer,parameter :: NCOMP = 3
    real :: start, finish
    double precision, dimension(GF%NGLLX) :: hxi, hpxi
    double precision, dimension(GF%NGLLY) :: heta, hpeta
    double precision, dimension(GF%NGLLZ) :: hgamma, hpgamma
    double precision :: h_x, h_y, h_z, h_xi, h_eta, h_gamma

    ! Initialize seismograms
    seismograms(:,:,:) = 0.d0

    if (DEBUG) call cpu_time(start)

    if (source%force) then
      call throwerror(-1, "Force source interpolation not implemented.")
    else

      ! Get interpolation values
      call lagrange_any(source%xi,GF%NGLLX,GF%xigll,hxi,hpxi)
      call lagrange_any(source%eta,GF%NGLLY,GF%yigll,heta,hpeta)
      call lagrange_any(source%gamma,GF%NGLLZ,GF%zigll,hgamma,hpgamma)

      do k=1, GF%NGLLZ
        do j=1, GF%NGLLY
          do i=1, GF%NGLLX

            iglob = GF%ibool(i,j,k, source%ispec)

            ! Compute derivative with respect to element coordinates
            h_xi = hpxi(i)* heta(j) * hgamma(k)
            h_eta = hxi(i) * hpeta(j) * hgamma(k)
            h_gamma = hxi(i) * heta(j) * hpgamma(k)

            ! Transform derivative using the Jacobian values
            h_x = h_xi * source%xix + h_eta * source%etax + h_gamma * source%gammax
            h_y = h_xi * source%xiy + h_eta * source%etay + h_gamma * source%gammay
            h_z = h_xi * source%xiz + h_eta * source%etaz + h_gamma * source%gammaz

            ! Compute everything at once
            seismograms(:,:,:) = seismograms(:,:,:) &
              + (GF%displacement(:,:,1,iglob,:) * h_x) * source%Mxx &
              + (GF%displacement(:,:,2,iglob,:) * h_y) * source%Myy &
              + (GF%displacement(:,:,3,iglob,:) * h_z) * source%Mzz &
              + (GF%displacement(:,:,2,iglob,:) * h_x &
                + GF%displacement(:,:,1,iglob,:) * h_y) * source%Mxy &
              + (GF%displacement(:,:,3,iglob,:) * h_x &
                + GF%displacement(:,:,1,iglob,:) * h_z) * source%Mxz &
              + (GF%displacement(:,:,3,iglob,:) * h_y &
                + GF%displacement(:,:,2,iglob,:) * h_z) * source%Myz

          enddo
        enddo
      enddo
    endif

    if (DEBUG) call cpu_time(finish)
    if (DEBUG) print '("interpolate_source took ",f6.3," seconds.")', finish-start

  end subroutine interpolate_source

end submodule