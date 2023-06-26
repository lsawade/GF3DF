program TestInterpolation
!
!  This program tests interpolation a timeseries, which here is
!  simply a sine function. The user can choose between linear and spline.
!
!

use gf3d, only: spline1d, interp1d, throwerror, get_args, printhello

implicit none


! Variables
double precision, allocatable :: x(:), y(:), h, xq(:), yq(:), yp
integer :: i        ! dummy iterator
integer :: Nd  ! number of data points
integer :: Nq = 100 ! number of inteerpolation points
integer :: choice   ! Set interpolation
character(len=20), dimension(:), allocatable :: args

! Get command line argument
args = get_args()

call printhello()

! Evaluate Args
if ((size(args) .ne. 2)) then
  write (*,*) " "
  write (*,*) " "
  write (*,*) "Usage:"
  write (*,*) "----------------------------"
  write (*,*) " "
  write (*,*) "  test-interp <Nd> <{1,2}>"
  write (*,*) " "
  write (*,*) " Nd -- Number of datapoints on [0,10]"
  write (*,*) " "
  write (*,*) "  1 -- Linear interpolation"
  write (*,*) "  2 -- Cubic spline interpolation"
  write (*,*) " "
  write (*,*) " "
  call throwerror(-1, " test-interp takes one argument.")
endif

! NUmber of datapoints
read(args(1),*) Nd

allocate(x(Nd), y(Nd))
x(:) = 0.d0 + ((/(i, i=1, Nd, 1)/)-1) * 10.d0 / Nd
y(:) = sin(x) ** 2

! print *, 'Cubic Spline Interpolation Demo'
allocate(xq(Nq), yq(Nq))

! Get interpolation vector
xq = x(1) + ((/(i, i=1, Nq, 1)/)-1) * ((x(size(x))-x(1))/(Nq-1))

! Do either
if (args(2) == "1") then
    call interp1d(x,y,xq,yq)
else if (args(2) == "2") then
    call spline1d(x,y,xq,yq)
else
    write(*,*) "Argument given: ", args(1)
    stop "Wrong argument. Choose 1 (=linear) or 2 (=spline)."
endif

! Print the interpolateed data to file
do i=1,Nq
    print '(g18.6,1x,g18.6)', xq(i), yq(i)
enddo

end program TestInterpolation