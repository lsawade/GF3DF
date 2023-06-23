module gf_types

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

end module gf_types