module gf3d_get_seismograms
  implicit none
  private
  public :: get_seismograms
contains
  function get_seismograms(GF, sources) result(superseismograms)

    use gf, only: t_GF
    use sources, only: t_source
    use stf, only: get_stf, stf_convolution, correct_hdur
    use interpolation, only: interpolateMT
    use source_location, only: locate_sources

    ! In
    type(t_GF), intent(in) :: GF
    type(t_source), dimension(:) :: sources

    ! Local
    integer :: isource, i, j, k, iglob
    double precision :: t0, tc, hdur_diff
    double precision, dimension(:), allocatable :: t, stf
    double precision, dimension(:,:,:), allocatable :: seismograms
    double precision, dimension(:,:,:), allocatable :: convolution
    double precision, dimension(:,:,:,:,:,:,:), allocatable :: displacement

    ! Out
    double precision, dimension(:,:,:), allocatable :: superseismograms

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

      ! Get displacement array
      do k=1,GF%ngllz
        do j=1,GF%nglly
          do i=1,GF%ngllx
            iglob = GF%ibool(i,j,k, sources(isource)%ispec)
            displacement(:,:,:,i,j,k, :) = GF%displacement(:,:,:,iglob,:)
          enddo
        enddo
      enddo

      ! Interpolate displacement array to seismograms
      call interpolateMT(&
        displacement, size(GF%displacement,1), &
        GF%nsteps, GF%ngllx, GF%nglly, GF%ngllz, GF%xigll, GF%yigll, GF%zigll, &
        sources(isource)%Mxx, sources(isource)%Myy, sources(isource)%Mzz, &
        sources(isource)%Mxy, sources(isource)%Mxz, sources(isource)%Myz, &
        sources(isource)%xi, sources(isource)%eta,sources(isource)%gamma, &
        sources(isource)%xix, sources(isource)%xiy, sources(isource)%xiz, &
        sources(isource)%etax, sources(isource)%etay, sources(isource)%etaz, &
        sources(isource)%gammax, sources(isource)%gammay, sources(isource)%gammaz, &
        seismograms)

      ! Convolve seismogram array with STF
      convolution(:,:,:) = stf_convolution(seismograms, stf, GF%dt, tc)

      superseismograms(:,:,:) = superseismograms(:,:,:) + convolution(:,:,:)

      write (*,*) 'Max disp: ', maxval(displacement)
      write (*,*) 'Max seis: ', maxval(seismograms)
      write (*,*) 'Max stf: ', maxval(stf)
      write (*,*) 'Max conv: ', maxval(convolution)
      write (*,*) 'Max supe: ', maxval(superseismograms)
    enddo


  end function get_seismograms


end module gf3d_get_seismograms