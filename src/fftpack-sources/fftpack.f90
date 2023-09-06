module fftpack
    use fftpack_kind

    implicit none
    private

    public :: zffti, zfftf, zfftb
    public :: fft, ifft
    public :: fftshift, ifftshift
    public :: fftfreq, rfftfreq

    public :: dffti, dfftf, dfftb
    public :: rfft, irfft

    public :: dzffti, dzfftf, dzfftb

    public :: dcosqi, dcosqf, dcosqb
    public :: dcosti, dcost
    public :: dct, idct
    public :: dct_t1i, dct_t1
    public :: dct_t23i, dct_t2, dct_t3

    public :: rk

    interface

        !> Version: experimental
        !>
        !> Initialize `zfftf` and `zfftb`.
        !> ([Specification](../page/specs/fftpack.html#zffti))
        pure subroutine zffti(n, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: wsave(*)
        end subroutine zffti

        !> Version: experimental
        !>
        !> Forward transform of a complex periodic sequence.
        !> ([Specification](../page/specs/fftpack.html#zfftf))
        pure subroutine zfftf(n, c, wsave)
            import rk
            integer, intent(in) :: n
            complex(kind=rk), intent(inout) :: c(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine zfftf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `zfftf`.
        !> ([Specification](../page/specs/fftpack.html#zfftb))
        pure subroutine zfftb(n, c, wsave)
            import rk
            integer, intent(in) :: n
            complex(kind=rk), intent(inout) :: c(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine zfftb

        !> Version: experimental
        !>
        !> Initialize `dfftf` and `dfftb`.
        !> ([Specification](../page/specs/fftpack.html#dffti))
        pure subroutine dffti(n, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: wsave(*)
        end subroutine dffti

        !> Version: experimental
        !>
        !> Forward transform of a real periodic sequence.
        !> ([Specification](../page/specs/fftpack.html#dfftf))
        pure subroutine dfftf(n, r, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(inout) :: r(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dfftf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `dfftf`.
        !> ([Specification](../page/specs/fftpack.html#dfftb))
        pure subroutine dfftb(n, r, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(inout) :: r(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dfftb

        !> Version: experimental
        !>
        !> Initialize `dzfftf` and `dzfftb`.
        !> ([Specification](../page/specs/fftpack.html#dzffti))
        pure subroutine dzffti(n, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: wsave(*)
        end subroutine dzffti

        !> Version: experimental
        !>
        !> Simplified forward transform of a real periodic sequence.
        !> ([Specification](../page/specs/fftpack.html#dzfftf))
        pure subroutine dzfftf(n, r, azero, a, b, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(in) :: r(*)
            real(kind=rk), intent(out) :: azero
            real(kind=rk), intent(out) :: a(*), b(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dzfftf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `dzfftf`.
        !> ([Specification](../page/specs/fftpack.html#dzfftb))
        pure subroutine dzfftb(n, r, azero, a, b, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: r(*)
            real(kind=rk), intent(in) :: azero
            real(kind=rk), intent(in) :: a(*), b(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dzfftb

        !> Version: experimental
        !>
        !> Initialize `dcosqf` and `dcosqb`.
        !> ([Specification](../page/specs/fftpack.html#initialize-dct-2-3-dcosqi-or-dct_t23i))
        pure subroutine dcosqi(n, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: wsave(*)
        end subroutine dcosqi

        !> Version: experimental
        !>
        !> Forward transform of quarter wave data.
        !> ([Specification](../page/specs/fftpack.html#compute-dct-3-dcosqf-or-dct_t3))
        pure subroutine dcosqf(n, x, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(inout) :: x(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dcosqf

        !> Version: experimental
        !>
        !> Unnormalized inverse of `dcosqf`.
        !> ([Specification](../page/specs/fftpack.html#compute-dct-2-dcosqb-or-dct_t2))
        pure subroutine dcosqb(n, x, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(inout) :: x(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dcosqb

        !> Version: experimental
        !>
        !> Initialize `dcost`.
        !> ([Specification](../page/specs/fftpack.html#initialize-dct-1-dcosti-or-dct_t1i))
        pure subroutine dcosti(n, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(out) :: wsave(*)
        end subroutine dcosti

        !> Version: experimental
        !>
        !> Discrete fourier cosine transform of an even sequence.
        !> ([Specification](../page/specs/fftpack.html#compute-dct-1-dcost-or-dct_t1))
        pure subroutine dcost(n, x, wsave)
            import rk
            integer, intent(in) :: n
            real(kind=rk), intent(inout) :: x(*)
            real(kind=rk), intent(in) :: wsave(*)
        end subroutine dcost

    end interface

    !> Version: experimental
    !>
    !> Forward transform of a complex periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#fft))
    interface fft
        module procedure fft_rk
    end interface fft

    !> Version: experimental
    !>
    !> Backward transform of a complex periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#ifft))
    interface ifft
        module procedure ifft_rk
    end interface ifft

    !> Version: experimental
    !>
    !> Forward transform of a real periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#rfft))
    interface rfft
        module procedure rfft_rk
    end interface rfft

    !> Version: experimental
    !>
    !> Backward transform of a real periodic sequence.
    !> ([Specifiction](../page/specs/fftpack.html#irfft))
    interface irfft
        module procedure irfft_rk
    end interface irfft

    !> Version: experimental
    !>
    !> Dsicrete cosine transforms.
    !> ([Specification](../page/specs/fftpack.html#simplified-dct-of-types-1-2-3-dct))
    interface dct
        module procedure dct_rk
    end interface dct

    !> Version: experimental
    !>
    !> Inverse discrete cosine transforms.
    !> ([Specification](../page/specs/fftpack.html#simplified-inverse-dct-of-types-1-2-3-idct))
    interface idct
        module procedure idct_rk
    end interface idct

    !> Version: experimental
    !>
    !> Initialize DCT type-1
    !> ([Specification](../page/specs/fftpack.html#initialize-dct-1-dcosti-or-dct_t1i))
    interface dct_t1i
        procedure :: dcosti
    end interface dct_t1i

    !> Version: experimental
    !>
    !> Perform DCT type-1
    !> ([Specification](../page/specs/fftpack.html#compute-dct-1-dcost-or-dct_t1))
    interface dct_t1
        procedure :: dcost
    end interface dct_t1

    !> Version: experimental
    !>
    !> Initialize DCT types 2, 3
    !> ([Specification](../page/specs/fftpack.html#initialize-dct-2-3-dcosqi-or-dct_t23i))
    interface dct_t23i
        procedure :: dcosqi
    end interface dct_t23i

    !> Version: experimental
    !>
    !> Perform DCT type-2
    !> ([Specification](../page/specs/fftpack.html#compute-dct-2-dcosqb-or-dct_t2))
    interface dct_t2
        procedure :: dcosqb
    end interface dct_t2

    !> Version: experimental
    !>
    !> Perform DCT type-3
    !> ([Specification](../page/specs/fftpack.html#compute-dct-3-dcosqf-or-dct_t3))
    interface dct_t3
        procedure :: dcosqf
    end interface dct_t3

    !> Version: experimental
    !>
    !> Shifts zero-frequency component to center of spectrum.
    !> ([Specifiction](../page/specs/fftpack.html#fftshift))
    interface fftshift
        module procedure fftshift_crk, fftshift_rrk
    end interface fftshift

    !> Version: experimental
    !>
    !> Shifts zero-frequency component to beginning of spectrum.
    !> ([Specifiction](../page/specs/fftpack.html#ifftshift))
    interface ifftshift
       module procedure ifftshift_crk, ifftshift_rrk
    end interface ifftshift


contains

!
! ----------------------------------------------------------------------------
!
    !> Shifts zero-frequency component to center of spectrum for `complex` type.
    pure function fftshift_crk(x) result(result)

        complex(kind=rk), intent(in) :: x(:)
        complex(kind=rk), dimension(size(x)) :: result

        result = cshift(x, shift=-floor(0.5_rk*size(x)))

    end function fftshift_crk

    !> Shifts zero-frequency component to center of spectrum for `real` type.
    pure function fftshift_rrk(x) result(result)

        real(kind=rk), intent(in) :: x(:)
        real(kind=rk), dimension(size(x)) :: result

        result = cshift(x, shift=-floor(0.5_rk*size(x)))

    end function fftshift_rrk

!
! ----------------------------------------------------------------------------
!


!> Shifts zero-frequency component to beginning of spectrum for `complex` type.
    pure function ifftshift_crk(x) result(result)
        complex(kind=rk), intent(in) :: x(:)
        complex(kind=rk), dimension(size(x)) :: result

        result = cshift(x, shift=-ceiling(0.5_rk*size(x)))

    end function ifftshift_crk

    !> Shifts zero-frequency component to beginning of spectrum for `real` type.
    pure module function ifftshift_rrk(x) result(result)
        real(kind=rk), intent(in) :: x(:)
        real(kind=rk), dimension(size(x)) :: result

        result = cshift(x, shift=-ceiling(0.5_rk*size(x)))

    end function ifftshift_rrk

!
! ----------------------------------------------------------------------------
!

 !> Discrete cosine transforms of types 1, 2, 3.
    pure function dct_rk(x, n, type) result(result)
        real(kind=rk), intent(in) :: x(:)
        integer, intent(in), optional :: n
        integer, intent(in), optional :: type
        real(kind=rk), allocatable :: result(:)

        integer :: lenseq, lensav, i
        real(kind=rk), allocatable :: wsave(:)

        if (present(n)) then
            lenseq = n
            if (lenseq <= size(x)) then
                result = x(:lenseq)
            else if (lenseq > size(x)) then
                result = [x, (0.0_rk, i=1, lenseq - size(x))]
            end if
        else
            lenseq = size(x)
            result = x
        end if

        ! Default to DCT-2
        if (.not.present(type)) then
            lensav = 3*lenseq + 15
            allocate (wsave(lensav))
            call dcosqi(lenseq, wsave)
            call dcosqb(lenseq, result, wsave)
            return
        end if

        if (type == 1) then  ! DCT-1
            lensav = 3*lenseq + 15
            allocate (wsave(lensav))
            call dcosti(lenseq, wsave)
            call dcost(lenseq, result, wsave)
        else if (type == 2) then  ! DCT-2
            lensav = 3*lenseq + 15
            allocate (wsave(lensav))
            call dcosqi(lenseq, wsave)
            call dcosqb(lenseq, result, wsave)
        else if (type == 3) then  ! DCT-3
            lensav = 3*lenseq + 15
            allocate (wsave(lensav))
            call dcosqi(lenseq, wsave)
            call dcosqf(lenseq, result, wsave)
        end if
    end function dct_rk

    !> Inverse discrete cosine transforms of types 1, 2, 3.
    pure function idct_rk(x, n, type) result(result)
        real(kind=rk), intent(in) :: x(:)
        integer, intent(in), optional :: n
        integer, intent(in), optional :: type
        real(kind=rk), allocatable :: result(:)

        integer :: lenseq, lensav, i
        real(kind=rk), allocatable :: wsave(:)

        if (present(n)) then
            lenseq = n
            if (lenseq <= size(x)) then
                result = x(:lenseq)
            else if (lenseq > size(x)) then
                result = [x, (0.0_rk, i=1, lenseq - size(x))]
            end if
        else
            lenseq = size(x)
            result = x
        end if

        ! Default to t=2; inverse DCT-2 is DCT-3
        if (.not.present(type)) then
            lensav = 3*lenseq + 15
            allocate (wsave(lensav))
            call dcosqi(lenseq, wsave)
            call dcosqf(lenseq, result, wsave)
            return
        end if

        if (type == 1) then  ! inverse DCT-1 is DCT-1
            lensav = 3*lenseq + 15
            allocate (wsave(lensav))
            call dcosti(lenseq, wsave)
            call dcost(lenseq, result, wsave)
        else if (type == 2) then  ! inverse DCT-2 is DCT-3
            lensav = 3*lenseq + 15
            allocate (wsave(lensav))
            call dcosqi(lenseq, wsave)
            call dcosqf(lenseq, result, wsave)
        else if (type == 3) then  ! inverse DCT-3 is DCT-2
            lensav = 3*lenseq + 15
            allocate (wsave(lensav))
            call dcosqi(lenseq, wsave)
            call dcosqb(lenseq, result, wsave)
        end if
    end function idct_rk

!
! ----------------------------------------------------------------------------
!

 !> Forward transform of a complex periodic sequence.
    pure function fft_rk(x, n) result(result)
        complex(kind=rk), intent(in) :: x(:)
        integer, intent(in), optional :: n
        complex(kind=rk), allocatable :: result(:)

        integer :: lenseq, lensav, i
        real(kind=rk), allocatable :: wsave(:)

        if (present(n)) then
            lenseq = n
            if (lenseq <= size(x)) then
                result = x(:lenseq)
            else if (lenseq > size(x)) then
                result = [x, ((0.0_rk, 0.0_rk), i=1, lenseq - size(x))]
            end if
        else
            lenseq = size(x)
            result = x
        end if

        !> Initialize FFT
        lensav = 4*lenseq + 15
        allocate (wsave(lensav))
        call zffti(lenseq, wsave)

        !> Forward transformation
        call zfftf(lenseq, result, wsave)

    end function fft_rk

!
! ----------------------------------------------------------------------------
!

    !> Backward transform of a complex periodic sequence.
    pure  function ifft_rk(x, n) result(result)
        complex(kind=rk), intent(in) :: x(:)
        integer, intent(in), optional :: n
        complex(kind=rk), allocatable :: result(:)

        integer :: lenseq, lensav, i
        real(kind=rk), allocatable :: wsave(:)

        if (present(n)) then
            lenseq = n
            if (lenseq <= size(x)) then
                result = x(:lenseq)
            else if (lenseq > size(x)) then
                result = [x, ((0.0_rk, 0.0_rk), i=1, lenseq - size(x))]
            end if
        else
            lenseq = size(x)
            result = x
        end if

        !> Initialize FFT
        lensav = 4*lenseq + 15
        allocate (wsave(lensav))
        call zffti(lenseq, wsave)

        !> Backward transformation
        call zfftb(lenseq, result, wsave)

    end function ifft_rk

!
! ----------------------------------------------------------------------------
!

    !> Forward transform of a real periodic sequence.
    pure module function rfft_rk(x, n) result(result)
        real(kind=rk), intent(in) :: x(:)
        integer, intent(in), optional :: n
        real(kind=rk), allocatable :: result(:)

        integer :: lenseq, lensav, i
        real(kind=rk), allocatable :: wsave(:)

        if (present(n)) then
            lenseq = n
            if (lenseq <= size(x)) then
                result = x(:lenseq)
            else if (lenseq > size(x)) then
                result = [x, (0.0_rk, i=1, lenseq - size(x))]
            end if
        else
            lenseq = size(x)
            result = x
        end if

        !> Initialize FFT
        lensav = 2*lenseq + 15
        allocate (wsave(lensav))
        call dffti(lenseq, wsave)

        !> Forward transformation
        call dfftf(lenseq, result, wsave)

    end function rfft_rk

!
! ----------------------------------------------------------------------------
!

    !> Backward transform of a real periodic sequence.
    pure module function irfft_rk(x, n) result(result)
        real(kind=rk), intent(in) :: x(:)
        integer, intent(in), optional :: n
        real(kind=rk), allocatable :: result(:)

        integer :: lenseq, lensav, i
        real(kind=rk), allocatable :: wsave(:)

        if (present(n)) then
            lenseq = n
            if (lenseq <= size(x)) then
                result = x(:lenseq)
            else if (lenseq > size(x)) then
                result = [x, (0.0_rk, i=1, lenseq - size(x))]
            end if
        else
            lenseq = size(x)
            result = x
        end if

        !> Initialize FFT
        lensav = 2*lenseq + 15
        allocate (wsave(lensav))
        call dffti(lenseq, wsave)

        !> Backward transformation
        call dfftb(lenseq, result, wsave)

    end function irfft_rk

!
! ----------------------------------------------------------------------------
!

    pure function fftfreq(n) result(out)
        integer, intent(in) :: n
        integer, dimension(n) :: out
        integer :: i

        out(1) = 0
        if (n == 1) return

        if (mod(n, 2) == 0) then  ! n even, smallest n = 2
            do i = 2, n/2
                out(i) = i-1
            end do
            out(n/2+1) = -n/2
            do i = n/2+2, n  ! only enters if n/2+2 <= n
                out(i) = out(i-1) + 1
            end do
        else  ! n odd, smallest n = 3
            do i = 2, n/2+1
                out(i) = i-1
            end do
            out(n/2+2) = -out(n/2+1)
            do i = n/2+3, n  ! only enters if n/2+3 <= n
                out(i) = out(i-1) + 1
            end do
        end if
    end function fftfreq

    !> Returns an integer array with the frequency values involved in the
    !> performed real FFT, ordered in the standard way (zero first, then
    !> positive frequencies, then, if applicable, the negative one).
    pure function rfftfreq(n) result(out)
        integer, intent(in) :: n
        integer, dimension(n) :: out
        integer :: i

        out(1) = 0
        if (n == 1) return

        if (mod(n,2) == 0) then  ! n even, smallest n = 2
            do i = 2, n-2, 2
                out(i) = out(i-1) + 1
                out(i+1) = out(i)
            end do
            out(n) = -n/2
        else  ! n odd, smallest n = 3
            do i = 2, n-1, 2
                out(i) = out(i-1) + 1
                out(i+1) = out(i)
            end do
        end if
    end function rfftfreq


end module fftpack
