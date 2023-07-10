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
    double precision, dimension(:,:,:,:,:,:,:) :: displacement
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

    ! In
    type(t_GF), intent(in) :: GF
    type(t_source), intent(in) :: source
    ! Out
    double precision, dimension(:,:,:) :: seismograms

    ! Local
    integer :: i,j,k,iglob
    real :: start, finish
    double precision, dimension(size(GF%displacement,1), 3, 3, GF%ngllx,GF%nglly,GF%ngllz,GF%nsteps) :: displacement

    if (DEBUG) call cpu_time(start)

    ! Get displacement array
    do k=1,GF%ngllz
      do j=1,GF%nglly
        do i=1,GF%ngllx

          displacement(:,:,:,i,j,k, :) =
        enddo
      enddo
    enddo

    if (source%force) then
      call throwerror(-1, "Force source interpolation not implemented.")
    else


      ! Get interpolation values
    call lagrange_any(xi,NGLLX,xigll,hxi,hpxi)
    call lagrange_any(eta,NGLLY,yigll,heta,hpeta)
    call lagrange_any(gamma,NGLLZ,zigll,hgamma,hpgamma)


    if (DEBUG) call cpu_time(startsub)

    seismograms(:,:,:) = 0.d0
    do k=1, NGLLZ
      do j=1, NGLLY
        do i=1, NGLLX

          iglob = GF%ibool(i,j,k, source%ispec)

          ! Compute derivative with respect to element coordinates
          h_xi = hpxi(i)* heta(j) * hgamma(k)
          h_eta = hxi(i) * hpeta(j) * hgamma(k)
          h_gamma = hxi(i) * heta(j) * hpgamma(k)

          ! Transform derivative using the Jacobian values
          h_x = h_xi * xix + h_eta * etax + h_gamma * gammax
          h_y = h_xi * xiy + h_eta * etay + h_gamma * gammay
          h_z = h_xi * xiz + h_eta * etaz + h_gamma * gammaz

          ! Compute everything at once
          seismograms(:,:,:) = seismograms(:,:,:) + (GF%displacement(:,:,1,iglob,:) * h_x) * moment_tensor(1) &
                                                  + (GF%displacement(:,:,2,iglob,:) * h_y) * moment_tensor(2) &
                                                  + (GF%displacement(:,:,3,iglob,:) * h_z) * moment_tensor(3) &
                                                  + 0.5 * (GF%displacement(:,:,2,iglob,:) * h_x &
                                                    + GF%displacement(:,:,1,iglob,:) * h_y) * moment_tensor(4) &
                                                  + 0.5 * (GF%displacement(:,:,3,iglob,:) * h_x &
                                                    + GF%displacement(:,:,1,iglob,:) * h_z) * moment_tensor(5) &
                                                  + 0.5 * (GF%displacement(:,:,3,iglob,:) * h_y &
                                                    + GF%displacement(:,:,2,iglob,:) * h_z) * moment_tensor(6) &
        enddo
      enddo
    enddo

    if (DEBUG) call cpu_time(finish)
    if (DEBUG) print '("interpolate_source took ",f6.3," seconds.")', finish-start

  end subroutine interpolate_source

end submodule