program write_source_time_functions
!
!  Auto source file test
!

  use gf3d, only: write_output_SAC
  use gf3d, only: get_stf, nextpower2
  use fftpack, only: fft, ifft, rk

  ! variable name
  logical :: OUTPUT_SEISMOS_SAC_ALPHANUM = .false.
  logical :: OUTPUT_SEISMOS_SAC_BINARY = .true.

  character(len=256) :: sisname
  character(len=256) :: OUTPUT_DIR = 'OUTPUT'
  character(len=256) :: model = "GLAD-M25"
  integer(kind=8), parameter :: NT = 4000
  integer(kind=8) :: i
  integer :: NP2
  double precision :: tc, dt, t0, hdur
  double precision, dimension(NT) :: t, stf
  complex(kind=rk), dimension(:), allocatable :: cstf
  double precision, dimension(:), allocatable :: acstf


  model = 'Test-Model'
  dt = 4.0
  hdur = 50.0
  tc = 200.0
  t0 = 0.0
  t(:) = t0 + ((/(i, i=1, NT, 1)/)-1) * dt

  call get_stf(t, tc, hdur, int(NT, kind=4), stf)

  write(sisname,"(a,'.',a)")
  write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
        "STF", "ERF", "TIM"

  call write_output_SAC(&
    stf, &
    1, &
    sisname, &
    "TIM", &
    1999, &
    365, &
    23, &
    59, &
    59.d0, &
    0.1d0, &
    tc, &
    "STF_TEST", &
    44.9d0, &
    99.9d0, &
    99.9d0, &
    9.9d0, &
    "ERF", &
    "STF", &
    44.9d0, &
    99.9d0, &
    0.0d0, &
    0.0d0, &
    dt, &
    0.0d0, &
    NT, &
    OUTPUT_SEISMOS_SAC_ALPHANUM, &
    OUTPUT_SEISMOS_SAC_BINARY, &
    model, &
    OUTPUT_DIR)

  NP2 = nextpower2(2 * int(NT, kind=4))


  cstf = fft(cmplx(stf,kind=rk), NP2)

  write(*,*) 'Sizes', NT, NP2

  allocate(acstf(NP2))
  acstf = dble(abs(cstf))

  write(sisname,"(a,'.',a)")
  write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
    "STF", "ERF", "FRE"

  call write_output_SAC(&
    stf, &
    1, &
    sisname, &
    "FRE", &
    1999, &
    365, &
    23, &
    59, &
    59.d0, &
    0.1d0, &
    tc, &
    "STF_TEST", &
    44.9d0, &
    99.9d0, &
    99.9d0, &
    9.9d0, &
    "ERF", &
    "STF", &
    44.9d0, &
    99.9d0, &
    0.0d0, &
    0.0d0, &
    dt, &
    0.0d0, &
    NT, &
    OUTPUT_SEISMOS_SAC_ALPHANUM, &
    OUTPUT_SEISMOS_SAC_BINARY, &
    model, &
    OUTPUT_DIR)


  cstf = ifft(cstf, NP2)

  stf = dble(real(cstf, kind=8))

  write(sisname,"(a,'.',a)")
  write(sisname,"('/',a,'.',a,'.',a3,'.sem')") &
    "STF", "ERF", "IFR"

  call write_output_SAC(&
    stf, &
    1, &
    sisname, &
    "IFR", &
    1999, &
    365, &
    23, &
    59, &
    59.d0, &
    0.1d0, &
    tc, &
    "STF_TEST", &
    44.9d0, &
    99.9d0, &
    99.9d0, &
    9.9d0, &
    "ERF", &
    "STF", &
    44.9d0, &
    99.9d0, &
    0.0d0, &
    0.0d0, &
    dt, &
    0.0d0, &
    NT, &
    OUTPUT_SEISMOS_SAC_ALPHANUM, &
    OUTPUT_SEISMOS_SAC_BINARY, &
    model, &
    OUTPUT_DIR)

  ! call cpu_time(finish)
  ! print '("Time = ",f6.3," seconds.")',finish-start
  ! print '("AvgTime = ",f6.3," seconds.")',(finish-start)/niter


 end program write_source_time_functions
