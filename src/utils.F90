module utils

  implicit none
  private
  public :: throwerror, scaleM, nextpower2, is_digit, is_numeric, get_args, &
            file_exists, dir_exists, dir_writeable, delete_file, mkdir

contains

  ! --------------------------------------------------------------------------

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

  ! --------------------------------------------------------------------------

  double precision function scaleM(M)

    use constants, only: RHOAV, R_PLANET, PI, GRAV

    double precision, intent(in) :: M

    scaleM = M / (1.d7 * RHOAV * (R_PLANET**5) * PI*GRAV*RHOAV)

  end function scaleM

  ! --------------------------------------------------------------------------

  pure integer function nextpower2(x_in)
    integer, intent(in)::x_in
    real ::x
    x = x_in
    nextpower2 = 2**(floor(log(x)/log(2.0))+1)
  end function nextpower2

  ! --------------------------------------------------------------------------

  logical function is_numeric(char)

    ! returns .true. if input character is a number

    implicit none
    character(len=1), intent(in) :: char

    is_numeric = .false.

    if ( index('0123456789', char) /= 0) then
      is_numeric = .true.
    endif

  end function

  !---------------------------------------------------------------------------

  logical function is_digit(char)

    ! returns .true. if input character is a number or a '.'

    implicit none
    character(len=1), intent(in) :: char

    is_digit = .false.

    if ( index('0123456789.', char) /= 0) then
      is_digit = .true.
    endif

  end function

  !---------------------------------------------------------------------------

  function get_args() result(args)

    integer :: ix, num_args                             !  some declarations
    character(len=20), dimension(:), allocatable :: args

    ! filename = '../single_element.h5' ! refer to line 5 of hdf5.js
    num_args = command_argument_count()

    ! I've omitted checking the return status of the allocation
    allocate(args(num_args))

    ! Actually read each argument.
    do ix = 1, num_args
      call get_command_argument(ix,args(ix))
    end do

  end function

  !---------------------------------------------------------------------------

  logical function dir_exists(dirname)

    implicit none
    character(len=*) :: dirname

    inquire (file=trim(dirname), exist=dir_exists)

  end function dir_exists

  !---------------------------------------------------------------------------

  logical function file_exists(filename)

    implicit none
    character(len=*) :: filename

    inquire (file=trim(filename), exist=file_exists)

  end function file_exists

  logical function dir_empty(dirname)

    implicit none
    character(len=*) :: dirname

  end function dir_empty

  !---------------------------------------------------------------------------

  function dir_writeable(dirname) result(res)

    implicit none
    character(len=*), intent(in) :: dirname
    integer :: unitno = 1234
    integer :: stat
    logical :: res
    character(len=30) :: testfile = 'deleteme.txt'

    ! Test whether the directory exists
    open(newunit=unitno,file=trim(dirname)//'/'//trim(testfile),status='replace',iostat=stat)

    if (stat .ne. 0) then

      res = .false.

    ! if exists close unit, set writeable to true and delete temp file
    else
      close (unitno)
      res = .true.

      call delete_file(trim(dirname)//'/'//trim(testfile))
    endif

  end function dir_writeable

  !---------------------------------------------------------------------------

  subroutine delete_file(filename)

    character(len=*) :: filename
    integer :: stat

    open(unit=1234, iostat=stat, file=filename, status='old')
    if (stat == 0) close(1234, status='delete')

  end subroutine

  subroutine mkdir(dirname)
    character(len=*) :: dirname

    call system('mkdir ' // trim(dirname))

  end subroutine mkdir

  ! logical function file_exists(filename)

  !   use iso_c_binding
  !   character(len=*) :: filename

  !   integer(c_int) :: mode, exist, time
  !   integer(c_int) :: true = 1
  !   integer(c_int) :: false = 0
  !   character(len=256, kind=c_char) :: cfilename


  !   ! Convert filename
  !   cfilename = trim(filename)

  !   ! Get file info
  !   call file_info(trim(cfilename), mode, exist, time)

  !   ! Convert c type to fortran type
  !   if (exist == true) then
  !     file_exists = .true.
  !   elseif (exist == false) then
  !     file_exists = .false.
  !   else
  !     write(*,*) "Something went wrong in file_exists."
  !     stop
  !   endif

  ! end function file_exists

  ! logical function dir_exists(dirname)

  !   use iso_c_binding
  !   character(len=*) :: dirname
  !   ! integer :: value
  !   integer(c_int) :: exist = 0
  !   integer(c_int) :: true = 1
  !   integer(c_int) :: false = 0
  !   character(len=256, kind=c_char) :: cdirname

  !   cdirname = dirname

  !   call dir_exists_c(cdirname, exist)

  !   ! Convert c type to fortran type
  !   if (exist == true) then
  !     dir_exists = .true.
  !   elseif (exist == false) then
  !     dir_exists = .false.
  !   endif

  ! end function


  ! subroutine dir_exists_c(dirname, exists) bind(c,name="dir_exists_c")
  !   use iso_c_binding
  !   character(kind=c_char), intent(in) :: dirname(*)
  !   integer(c_int), intent(out) :: exists
  ! end subroutine

  ! subroutine file_info(filename,mode,exist,time) bind(c,name="file_info")
  !   use iso_c_binding
  !   character(kind=c_char),intent(in) :: filename(*)
  !   integer(c_int),intent(out) :: mode,exist,time
  ! end subroutine

end module