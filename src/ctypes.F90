module ctypes

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
      character(len=5), dimension(:),         allocatable  :: networks, stations
      real(kind=8),     dimension(:,:,:,:,:), allocatable  :: latitudes, longitudes

      ! Arrays
      real(kind=8),     dimension(:),         allocatable  :: rspl, ellipticity_spline, ellipticity_spline2
      integer,          dimension(:,:),       allocatable  :: bathy
      real(kind=8),     dimension(:,:,:,:,:), allocatable  :: displacement
      real(kind=8),     dimension(:,:),       allocatable  :: xyz
      integer(kind=8),  dimension(:,:,:,:),   allocatable  :: ibool
      integer(kind=8),  dimension(:),         allocatable  :: adjacency
      integer(kind=8),  dimension(:),         allocatable  :: xadj

      ! number of nodes for 2D and 3D shape functions for hexahedra with 27 nodes
      integer :: anchor_iax(27),anchor_iay(27),anchor_iaz(27)

    end type

  type, public :: t_MT

    ! Geographical Moment tensor
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

    ! Moment magnitude and scalar moment
    double precision :: Mw = 0.d0
    double precision :: M0 = 0.d0

    ! Local element coordinates
    double precision :: xi_source = 0.d0
    double precision :: eta_source = 0.d0
    double precision :: gamma_source = 0.d0
    double precision :: tshift_src = 0.d0
    double precision :: hdur = 0.d0
    double precision :: hdur_Gaussian = 0.d0

    ! Source coordinates in spherical coordinates
    double precision :: theta_source = 0.d0
    double precision :: phi_source = 0.d0
    double precision :: depth_source = 0.d0

  end type

  type, public :: t_F

    ! parameters for a force source located exactly at a grid point
    integer :: force_stf
    double precision :: factor_force_source
    double precision :: comp_dir_vect_source_E
    double precision :: comp_dir_vect_source_N
    double precision :: comp_dir_vect_source_Z_UP

    ! In-mesh, cartesian moment tensor
    double precision :: Mxx = 0.d0
    double precision :: Myy = 0.d0
    double precision :: Mzz = 0.d0
    double precision :: Mxy = 0.d0
    double precision :: Mxz = 0.d0
    double precision :: Myz = 0.d0

    ! Moment magnitude and scalar moment
    double precision :: Mw = 0.d0
    double precision :: M0 = 0.d0

    ! Local element coordinates
    double precision :: xi_source = 0.d0
    double precision :: eta_source = 0.d0
    double precision :: gamma_source = 0.d0
    double precision :: tshift_src = 0.d0
    double precision :: hdur = 0.d0
    double precision :: hdur_Gaussian = 0.d0

    ! Source coordinates in spherical coordinates
    double precision :: theta_source = 0.d0
    double precision :: phi_source = 0.d0
    double precision :: depth_source = 0.d0

  end type

  type, public :: t_Sources

    integer :: N_MT
    integer :: N_F

    type(t_MT), dimension(:), allocatable :: moment_tensors
    type(t_F), dimension(:), allocatable :: force_sources

  end type
end module