module sources

  use constants, only: HUGEVAL

  private
  public :: &
    t_source, &
    get_cmt_scalar_moment, &
    get_cmt_moment_magnitude, &
    get_cmt_moment_magnitude_from_M0, &
    read_cmt, &
    print_source


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

    ! PDE parameters
    character(len=30) :: pde_desc = " PDE "
    double precision :: pde_lat = 0.0
    double precision :: pde_lon = 0.0
    double precision :: pde_depth = 0.0
    double precision :: pde_mb = 0.0
    double precision :: pde_ms = 0.0
    character(len=30) :: pde_region = "------"

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


  interface
    double precision module function get_cmt_scalar_moment(Mxx,Myy,Mzz,Mxy,Mxz,Myz)
      double precision, intent(in) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
    end function get_cmt_scalar_moment

    double precision module function get_cmt_moment_magnitude(Mxx,Myy,Mzz,Mxy,Mxz,Myz)
      double precision, intent(in) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
    end function get_cmt_moment_magnitude

    double precision module function get_cmt_moment_magnitude_from_M0(M0)
      double precision, intent(in) :: M0
    end function

    module subroutine print_source(source, which)
      type(t_source) :: source
      integer :: which ! 1 cmtsolution, 2 meshinfo, 3 both
    end subroutine print_source

    module function read_cmt(filename) result(sources)
      character(len=*) :: filename
      type(t_source), dimension(:), allocatable :: sources
    end function

  end interface

end module sources