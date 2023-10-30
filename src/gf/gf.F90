module gf

  use constants, only: HUGEVAL
  use calendar, only: julian_day
  use hdf5, only: hid_t

  private
  public :: &
    read_GF, &
    print_GF, &
    t_GF, &
    free_GF

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
    real(kind=4),     dimension(:), allocatable  :: latitudes, longitudes
    real(kind=4),     dimension(:), allocatable  :: burials

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

    contains

      procedure :: get_kdtree => setup_point_search_arrays

  end type


  interface

    type(t_GF) module function read_GF(filename) result(GF)
      integer(hid_t)    :: file_id         ! file identifierd
      character(len=65) :: filename ! input variable
    end function

    module subroutine print_GF(GF)
      type(t_GF), intent(in) :: GF
    end subroutine

    module subroutine load_header(file_id, GF)
      integer(hid_t), intent(in)    :: file_id
      type(t_GF), intent(inout) :: GF
    end subroutine

    module subroutine load_arrays(file_id, GF)
      integer(hid_t), intent(in)    :: file_id
      type(t_GF), intent(inout) :: GF
    end subroutine

    module subroutine free_GF(GF)
      type(t_GF), intent(inout) :: GF
    end subroutine

  end interface





contains

  subroutine setup_point_search_arrays(GF)


    use constants, only: &
      IMAIN, NDIM

    use kdtree_search, only: kdtree_setup,kdtree_set_verbose, &
      kdtree_num_nodes,kdtree_nodes_location,kdtree_nodes_index

    implicit none

    ! variable types
    class(t_GF), intent(in)     :: GF

    ! local parameters
    integer(kind=8) :: iglob
    integer :: ispec,ier
    integer :: inodes
    ! determines tree points
    ! Pointers to variables, so that I don't have to change anything.

    ! define (i,j,k) indices of the control/anchor points
    call hex_nodes_anchor_ijk(&
        GF%ngllx,GF%nglly,GF%ngllz, &
        GF%midx,GF%midy,GF%midz, &
        GF%anchor_iax, GF%anchor_iay, GF%anchor_iaz)

    ! kd-tree setup for point localization
    !
    ! determines tree size
    kdtree_num_nodes = GF%nspec

    ! allocates tree arrays
    allocate(kdtree_nodes_location(NDIM,kdtree_num_nodes),stat=ier)
    if (ier /= 0) stop 'Error allocating kdtree_nodes_location arrays'
    allocate(kdtree_nodes_index(kdtree_num_nodes),stat=ier)
    if (ier /= 0) stop 'Error allocating kdtree_nodes_index arrays'

    ! prepares search arrays, each element takes its internal GLL points for tree search
    kdtree_nodes_index(:) = 0
    kdtree_nodes_location(:,:) = 0.0
    ! adds tree nodes
    inodes = 0

    ! sets up tree nodes
    do ispec = 1,GF%nspec
      iglob = GF%ibool(GF%midz, GF%midy, GF%midz,ispec)

      ! counts nodes
      inodes = inodes + 1
      if (inodes > kdtree_num_nodes) stop 'Error index inodes bigger than kdtree_num_nodes'

      ! adds node index (index points to same ispec for all internal GLL points)
      kdtree_nodes_index(inodes) = ispec

      ! adds node location
      kdtree_nodes_location(1,inodes) = GF%xyz(iglob, 1)
      kdtree_nodes_location(2,inodes) = GF%xyz(iglob, 2)
      kdtree_nodes_location(3,inodes) = GF%xyz(iglob, 3)
    enddo

    if (inodes /= kdtree_num_nodes) stop 'Error index inodes does not match kdtree_num_nodes'

    ! tree verbosity
    call kdtree_set_verbose(IMAIN)

    ! creates kd-tree for searching point locations in locate_point() routine
    call kdtree_setup()

    end subroutine setup_point_search_arrays


end module gf