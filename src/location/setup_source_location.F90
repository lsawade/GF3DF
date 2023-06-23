
module setup_source_location
  implicit none
  private
  public :: setup_point_search_arrays

contains

  subroutine setup_point_search_arrays(GF)

    use gf, only: t_GF

    use constants, only: &
      IMAIN, NDIM

    use kdtree_search, only: kdtree_setup,kdtree_set_verbose, &
      kdtree_num_nodes,kdtree_nodes_location,kdtree_nodes_index

    implicit none

    ! variable types
    type(t_GF), intent(in)     :: GF

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

end module