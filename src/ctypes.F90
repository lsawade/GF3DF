module ctypes

  use constants, only: HUGEVAL
  use calendar, only: julian_day
  implicit none
  private

  type, public :: t_GF

      ! Header parameter
      integer                        :: ellipticity = 0
      integer                        :: topography = 0
      integer                        :: do_adjacency_search = 0
      integer(kind=8)                :: nsteps = 0
      integer                        :: nspec = 0
      integer                        :: ngllx = 0
      integer                        :: nglly = 0
      integer                        :: ngllz = 0
      integer                        :: midx = 0
      integer                        :: midy = 0
      integer                        :: midz = 0
      integer                        :: nx_bathy = 0
      integer                        :: ny_bathy = 0
      real(kind=8)                   :: factor = 0
      real(kind=8)                   :: dt = 0
      real(kind=8)                   :: hdur = 0
      real(kind=8)                   :: tc = 0
      real(kind=8)                   :: resolution_topo_file = 0
      character(len=8), dimension(:), allocatable  :: networks
      character(len=32), dimension(:), allocatable  :: stations
      real(kind=8),     dimension(:), allocatable  :: latitudes, longitudes
      real(kind=8),     dimension(:), allocatable  :: burials

      ! Interpolation points, weights, values, and derivatives
      double precision, dimension(:), allocatable :: xigll, yigll, zigll
      double precision, dimension(:), allocatable :: wxgll, wygll, wzgll

      ! Arrays
      real(kind=8),     dimension(:),         allocatable  :: rspl, ellipticity_spline, ellipticity_spline2
      integer,          dimension(:,:),       allocatable  :: bathy
      real(kind=8),     dimension(:,:,:,:,:), allocatable  :: displacement
      real(kind=4),     dimension(:,:),       allocatable  :: xyz
      integer(kind=8),  dimension(:,:,:,:),   allocatable  :: ibool
      integer,  dimension(:),                 allocatable  :: adjacency
      integer,  dimension(:),                 allocatable  :: xadj

      ! number of nodes for 2D and 3D shape functions for hexahedra with 27 nodes
      integer :: anchor_iax(27),anchor_iay(27),anchor_iaz(27)

    end type

  type, public :: t_source

    ! Whether Force or moment tensor source
    logical :: force = .false.

    ! Eventname
    character(len=30) :: eventname

    ! Origin time parameters
    integer :: year = 1999
    integer :: month = 1
    integer :: day = 1
    integer :: jda = 1
    integer :: hour = 0
    integer :: minute = 0
    double precision :: second = 0.d0

    ! Geographical parameters
    double precision :: colatitude = 0.d0
    double precision :: latitude = 0.d0
    double precision :: longitude = 0.d0
    double precision :: depth = 0.d0
    double precision :: Mrr = 0.d0
    double precision :: Mtt = 0.d0
    double precision :: Mpp = 0.d0
    double precision :: Mrt = 0.d0
    double precision :: Mrp = 0.d0
    double precision :: Mtp = 0.d0

    ! In-mesh, cartesian moment tensor
    double precision :: Mxx = 0.d0
    double precision :: Myy = 0.d0
    double precision :: Mzz = 0.d0
    double precision :: Mxy = 0.d0
    double precision :: Mxz = 0.d0
    double precision :: Myz = 0.d0

    double precision, dimension(6) :: tensor = (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)

    ! Moment magnitude and scalar moment
    double precision :: Mw = 0.d0
    double precision :: M0 = 0.d0

    ! Mesh data
    integer :: ispec = 0
    double precision :: xi = 0.d0
    double precision :: eta = 0.d0
    double precision :: gamma = 0.d0
    double precision :: xix = 0.d0, xiy = 0.d0, xiz = 0.d0
    double precision :: etax = 0.d0, etay = 0.d0, etaz = 0.d0
    double precision :: gammax = 0.d0, gammay = 0.d0, gammaz = 0.d0

    ! STF parameters
    double precision :: time_shift = 0.d0
    double precision :: hdur = 0.d0
    double precision :: hdur_Gaussian = 0.d0

    ! Cartesian coordinates
    double precision :: x = 0.d0
    double precision :: y = 0.d0
    double precision :: z = 0.d0
    double precision :: x_found = 0.d0
    double precision :: y_found = 0.d0
    double precision :: z_found = 0.d0
    double precision :: x_target = 0.d0
    double precision :: y_target = 0.d0
    double precision :: z_target = 0.d0
    double precision :: final_distance = HUGEVAL

    double precision, dimension(3,3) :: nu


    ! Source coordinates in spherical coordinates
    double precision :: theta = 0.d0
    double precision :: phi = 0.d0
    double precision :: r_target = 0.d0
    double precision :: r_found = 0.d0
    double precision :: r0 = 0.d0
    double precision :: radius = 0.d0

    ! parameters for a force source located exactly at a grid point
    integer :: force_stf = 0
    double precision :: factor_force_source = 0.d0
    double precision :: comp_dir_vect_source_E = 0.d0
    double precision :: comp_dir_vect_source_N = 0.d0
    double precision :: comp_dir_vect_source_Z_UP = 0.d0

  ! contains
  !     private
  !     procedure, pass :: write => write_source
  !     ! procedure, pass :: read => read_source
  !     generic, public :: write(formatted) => write
  !     ! generic, public :: read(formatted)  => read

  end type

  contains

    subroutine write_source(source, unit, iotype, v_list, iostat, iomsg)

    class(t_source), intent(in) :: source
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in)  :: v_list(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg

    if (source%force) then
      ! N = f"{self.force_no:d}".zfill(3)
      ! rstr = f"FORCE  {N}\n"
      ! rstr += f"time shift:{self.time_shift:15.4f}    ! s\n"
      ! rstr += f"half duration:{self.hdur:9.1f}       ! Half duration (s) for Gaussian/Step function, frequency (Hz) for Ricker\n"
      ! rstr += f"latitude:{self.latitude:17.4f}    ! Degree \n"
      ! rstr += f"longitude:{self.longitude:16.4f}    ! Degree\n"
      ! rstr += f"depth:{self.depth:20.4f}    ! km\n"
      ! rstr += f"source time function:{self.stf:2d}       ! 0=Gaussian function, 1=Ricker wavelet, 2=Step function\n"
      ! FF = float_to_str(self.forcefactor, 5)
      ! E = float_to_str(self.vector_E, 3)
      ! N = float_to_str(self.vector_N, 3)
      ! Z = float_to_str(self.vector_Z_UP, 3)
      ! rstr += f"factor force source: {FF:<8} ! Newton\n"
      ! rstr += f"component dir vect source E: {E:>9}\n"
      ! rstr += f"component dir vect source N: {N:>9}\n"
      ! rstr += f"component dir vect source Z_UP: {Z:>5}  ! Upward force\n"
    else

      ! hello this is string for printing source.

    endif


    end subroutine write_source


end module