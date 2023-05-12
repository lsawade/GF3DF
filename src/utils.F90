module utils

  implicit none
  private
  public :: throwerror

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


end module