module utils

  implicit none
  private
  public :: throwerror, scaleM

contains

  subroutine throwerror(errorflag, message)

    ! In
    integer, intent(in) :: errorflag
    character(len=*), intent(in) :: message

    ! Local
    integer :: i

    if (errorflag /= 0) then
        write (*,*) message
        write (*,*) "Error code: ", errorflag

        ! i=1/(errorflag-errorflag)
        call exit()
    end if

  end subroutine throwerror

  double precision function scaleM(M)

    use constants, only: RHOAV, R_PLANET, PI, GRAV

    double precision, intent(in) :: M

    scaleM = M / (1.d7 * RHOAV * (R_PLANET**5) * PI*GRAV*RHOAV)

  end function scaleM

end module