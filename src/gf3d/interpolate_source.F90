submodule (gf3d:seismograms) source_interpolation

  implicit none

contains

  module subroutine interpolate_source(GF, source, seismograms)

    use gf, only: t_GF
    use sources, only: t_source
    use interpolation, only: interpolateMT

    ! In
    type(t_GF), intent(in) :: GF
    type(t_source), intent(in) :: source
    ! Out
    double precision, dimension(:,:,:) :: seismograms

    ! Local
    integer :: i,j,k,iglob
    double precision, dimension(:,:,:,:,:,:,:), allocatable :: displacement

    ! Allocate temporary displacement array
    allocate(displacement(size(GF%displacement,1), 3, 3, GF%ngllx,GF%nglly,GF%ngllz,GF%nsteps))

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
  end subroutine interpolate_source

end submodule