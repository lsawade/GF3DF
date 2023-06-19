module utils

  implicit none
  private
  public :: throwerror, scaleM, nextpower2

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

  pure integer  function nextpower2(x_in)
    integer, intent(in)::x_in
    real ::x
    x = x_in
    nextpower2 = 2**(floor(log(x)/log(2.0))+1)
  end function nextpower2

end module