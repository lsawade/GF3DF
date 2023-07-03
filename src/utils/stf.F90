

module stf
!
!  This module contains all functions related to source time function
!  generation, convolution and parameter evaluation.
!
    implicit none
    private
    public :: get_stf, stf_convolution, correct_hdur


contains

  subroutine get_stf(t, tc, hdur, NT, stf)
  !
  !  Computes STF using an Error function.
  !
  !  This function computes a step function to be convolved with Green function
  !  data. It uses an input time vector time-center and half-duration and
  !  populates an array with the output STF
  !
    use error_functions, only: netlib_specfun_erf

    double precision, dimension(NT), intent(in) :: t        ! Time vector
    double precision, intent(in) :: tc                      ! center of step in error function
    double precision, intent(in) :: hdur                    ! half duration
    integer :: NT                                           ! Length of timevector
    double precision, dimension(NT), intent(out) :: stf     ! Output source time function

    integer :: i

    do i=1,NT
       stf(i) = 0.5*netlib_specfun_erf((t(i)-tc)/hdur)+0.5
    enddo

  end subroutine get_stf

  function stf_convolution(data, stf, dt, tc)
  !
  !  Convolves a set of seismograms with a source time function.
  !
  !  The function takes an 3D array of data of which the last dimension are the
  !  time series to be convolved with an equal length STF. The convolved traces
  !  are shifted back in time to make the convolution independent of possible
  !  time shifts

    use fftpack, only: rk, fft, ifft, fftfreq
    use utils, only: nextpower2
    use constants, only: PI, DEBUG

    ! In
    double precision, dimension(:,:,:), intent(in) :: data             ! array
    double precision, dimension(:), intent(in) :: stf                  ! STF
    double precision, intent(in) :: dt                                 ! sample spacing
    double precision, intent(in) :: tc                                 ! time shift

    ! Local
    integer :: k, Nk, i, Ni, NP2, NT
    real :: start, finish
    complex(kind=rk), dimension(:), allocatable :: cstf
    complex(kind=rk), dimension(:), allocatable :: pshift
    complex(kind=rk), dimension(:), allocatable :: convo
    complex(kind=rk), dimension(:), allocatable :: cseis
    double precision, dimension(:), allocatable :: outseis
    double precision, dimension(:), allocatable :: freqs

    ! Out
    double precision, dimension(:,:,:), allocatable :: stf_convolution ! Output traces

    !------

    if (DEBUG) call cpu_time(start)


    ! Get shape of stf
    Nk = size(data, dim=1)
    Ni = size(data, dim=2)
    NT = size(data, dim=3)

    ! Get doubled next power of 2 (to be safe)
    NP2 = nextpower2(2 * int(NT, kind=4))

    ! Allocate
    allocate(freqs(NP2), pshift(NP2))
    allocate(convo(NP2))
    allocate(cseis(NT))
    allocate(outseis(NP2))
    allocate(stf_convolution(Nk,Ni,NT))

    ! Get complex STF
    cstf = fft(cmplx(stf,kind=rk), NP2)

    ! Get frequencies associated with the fourier transofrm to compute phase shift
    freqs = dble(fftfreq(NP2)) / (dt * NP2)

    ! Create phase shift to correct convolution timeshift
    pshift(:) = exp(complex(0,-1) * (-tc) * freqs * 2 * PI)

    ! Loop over traces to convolve
    do k=1,Nk
      do i=1,Ni

        ! Transfrom seismogram to complex
        cseis(:) = cmplx(data(k,i,:), kind=rk)

        ! Convolve seismogram STF and phaseshift
        convo(:) = fft(cseis, n=NP2) * dt * cstf * pshift

        ! IFFT
        outseis(:) = dble(ifft(convo, n=NP2)) / NP2

        ! Extract shorter original shape
        stf_convolution(k,i,:) = outseis(1:NT)

      enddo
    enddo


    if (DEBUG) call cpu_time(finish)
    if (DEBUG) print '("\n    stf_convolution took ",f6.3," seconds.\n")', finish-start


  end function stf_convolution


  function correct_hdur(hdur_cmt, hdur_gf)
  !
  !  Calculates corrected half duration for STF convolution
  !
  !  This function takes in the requested half duration from a CMTSOLUTION
  !  and the half duration used for computing the database Green functions, and
  !  outputs the correct half duration for the Source Time function to convolve
  !  the Green functions with. Should the hdur of the CMT be too short the
  !  function output a half duration close to 0, and suggests a half duration
  !  for a Gaussian to convolve the data with after convolution with STF.
  !

    ! In
    double precision, intent(in) :: hdur_cmt        ! Requested hdur from CMT
    double precision, intent(in) :: hdur_gf         ! hdur from database

    ! Local
    double precision :: hdur_conv                   ! suggested hdur if cmt hdur too short

    ! Out
    double precision :: correct_hdur   ! OUTPUT half duration


    if ((hdur_cmt / 1.628) ** 2 .le. hdur_gf**2) then

        correct_hdur = 0.000001

        write (*,*) &
            "Requested half duration smaller than what was simulated.\n", &
            "Half duration set to ", correct_hdur," s to simulate a Heaviside function."

        hdur_conv = sqrt(hdur_gf**2 - (hdur_cmt / 1.628)**2)

        write (*,*) "Try convolving your seismogram with a Gaussian with ", &
                    hdur_conv, " standard deviation."
    else

      correct_hdur = sqrt((hdur_cmt / 1.628)**2 - hdur_gf**2)

    endif

  end function correct_hdur

end module stf