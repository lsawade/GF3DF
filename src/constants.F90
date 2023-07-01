module constants

  logical :: DEBUG = .false.

  ! Some constant numbers
  double precision, parameter :: ONE = 1.0
  double precision, parameter :: ZERO = 0.0
  double precision, parameter :: HALF = 0.5
  double precision, parameter :: TWO = 2.0
  double precision, parameter :: TINYVAL = 1.0E-9
  double precision, parameter :: HUGEVAL = 1.0E30
  double precision, parameter :: VERYSMALLVAL = 1E-24

  ! Degrees and Radians and their conversions
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: TWO_PI = 2.d0 * PI
  double precision, parameter :: PI_OVER_FOUR = PI / 4.d0
  double precision, parameter :: PI_OVER_TWO = PI / 2.0d0

  ! to convert angles from degrees to radians
  double precision, parameter :: DEGREES_TO_RADIANS = PI / 180.d0

  ! to convert angles from radians to degrees
  double precision, parameter :: RADIANS_TO_DEGREES = 180.d0 / PI

  ! Radius of unit sphere
  double precision, parameter :: R_UNIT_SPHERE = ONE

  ! small tolerance for conversion from x y z to r theta phi
  double precision, parameter :: SMALL_VAL_ANGLE = 1.0E-10

  ! Some fractions
  double precision, parameter :: ONE_THIRD   = 1.d0/3.d0
  double precision, parameter :: TWO_THIRDS  = 2.d0/3.d0
  double precision, parameter :: FOUR_THIRDS = 4.d0/3.d0

  ! maximum length of strings used for paths, reading from files, etc.
  integer, parameter :: MAX_STRING_LEN = 512

  ! input, output and main MPI I/O files
  ! note: careful with these unit numbers, we mostly use units in the 40-50 range.
  !       Cray Fortran e.g. reserves 0,5,6 (standard error,input,output units) and 100,101,102 (input,output,error unit)
  integer, parameter :: ISTANDARD_OUTPUT = 6     ! or for cray: 101
  ! I/O unit for file input,output
  integer, parameter :: IIN = 40,IOUT = 41

  ! uncomment this to write messages to a text file
  ! integer, parameter :: IMAIN = 42
  character(len=256) :: LOGFILE = 'gf3d_log.txt'

  ! uncomment this to write messages to the screen (slows down the code)
  integer, parameter :: IMAIN = ISTANDARD_OUTPUT

  ! I/O unit for sac files
  integer, parameter :: IOUT_SAC = 48

  ! 3D simulation/mesh
  integer, parameter :: NDIM = 3

  ! for the Gauss-Lobatto-Legendre points and weights
  double precision, parameter :: GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0

  ! number of iterations to solve the system for xi and eta
  ! setting it to 5 instead of 4 ensures that the result obtained is not compiler dependent
  ! (when using 4 only some small discrepancies were observed)
  integer, parameter :: NUM_ITER = 5

  ! For in-element point location
  logical, parameter :: USE_DISTANCE_CRITERION = .false.

  ! Earth related constants
  double precision, parameter :: EARTH_FLATTENING_F = 1.0 / 299.80
  double precision, parameter :: EARTH_ONE_MINUS_F_SQUARED = (1.0 - EARTH_FLATTENING_F)**2

  ! EARTH_R is the radius of the bottom of the oceans(radius of Earth in m)
  double precision, parameter :: EARTH_R = 6371000.0

  ! and in kilometers:
  double precision, parameter :: EARTH_R_KM = EARTH_R / 1000.0

  ! average density in the full Earth to normalize equation
  double precision, parameter :: EARTH_RHOAV = 5514.3

  ! standard gravity at the surface of the Earth
  double precision, parameter :: EARTH_STANDARD_GRAVITY = 9.80665  ! in m.s-2

  ! Even though there is ellipticity, use spherical Earth assumption for the
  ! conversion from geographical to spherical coordinates.
  logical, parameter :: ASSUME_PERFECT_SPHERE = .false.

  ! gravitational constant in S.I. units i.e. in m3 kg-1 s-2, or equivalently in N.(m/kg)^2
  ! DK DK April 2014: switched to the 2010 Committee on Data for Science and Technology (CODATA) recommended value
  double precision, parameter :: GRAV = 6.67384e-11  ! CODATA 2010

  ! FOR NOW USE STATIC --> planets could be added later or specfem constants
  ! could be dumped.
  double precision, parameter :: ONE_MINUS_F_SQUARED = EARTH_ONE_MINUS_F_SQUARED
  double precision, parameter :: R_PLANET = EARTH_R
  double precision, parameter :: R_PLANET_KM = EARTH_R_KM
  double precision, parameter :: RHOAV = EARTH_RHOAV

  ! Node setup
  integer, parameter :: NGNOD  = 27

  ! Conversion for the source triangle
  double precision, parameter :: SOURCE_DECAY_MIMIC_TRIANGLE = 1.628

  ! Shift used to convolved STF with
  double precision, parameter :: STF_SHIFT = 200.d0

  ! Parameter used in coordinate conversion. Kept to not change
  ! the conversion codes too much
  logical, parameter :: USE_OLD_VERSION_5_1_5_FORMAT = .false.

  ! Output related constants
  integer :: NCHANNELS = 3
  character(len=4), dimension(3) :: channels = (/"MXN", "MXE", "MXZ"/)
  integer, dimension(3) :: orientation = (/1,2,3/)

  ! maximum length of station and network name for receivers
  integer, parameter :: MAX_LENGTH_STATION_NAME = 32
  integer, parameter :: MAX_LENGTH_NETWORK_NAME = 8

  ! Lines in the source files
  integer, parameter :: NLINES_PER_CMTSOLUTION_SOURCE = 13
  integer, parameter :: NLINES_PER_FORCESOLUTION_SOURCE = 11

  ! Perturbation values
  character(len=4), dimension(11) :: partialnames = &
    (/"mrr", "mtt", "mpp", "mrt", "mrp", "mtp", "lat", "lon", "dep", "cmt", "hdr"/)
  double precision,parameter :: dmom  = 1.0d23
  double precision,parameter :: dlat  = 0.0001d0
  double precision,parameter :: dlon  = 0.0001d0
  double precision,parameter :: ddep  = 0.01d0
  double precision,parameter :: dcmt  = -1.d0
  double precision,parameter :: dhdur = 0.001d0

end module constants