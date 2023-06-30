!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

module point_location

  private
  public :: locate_point

contains

  subroutine locate_point(&
        x_target,y_target,z_target, lat_target,lon_target, ispec_selected, &
        nspec, NGLLX, NGLLY,NGLLZ,MIDX,MIDY,MIDZ, &
        ibool,xstore,ystore,zstore,xadj,adjncy, &
        xigll, yigll, zigll, anchor_iax, anchor_iay, anchor_iaz, &
        POINT_CAN_BE_BURIED, &
        xi,eta,gamma, &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,x,y,z, &
        distmin_not_squared)

    ! locates target point inside best mesh element
    use constants, only: HUGEVAL, USE_DISTANCE_CRITERION, R_PLANET_KM, DEBUG, IMAIN

    ! for point search
    use kdtree_search, only: kdtree_find_nearest_neighbor

    implicit none

    ! IN
    double precision,intent(in) :: x_target,y_target,z_target
    double precision,intent(in) :: lat_target,lon_target

    integer, intent(in) :: nspec, NGLLX, NGLLY, NGLLZ, MIDX, MIDY, MIDZ
    integer(kind=8), dimension(:,:,:,:), intent(in) :: ibool
    integer, dimension(:), intent(in) :: xadj, adjncy
    real(kind=4), dimension(:), intent(in) :: xstore, ystore, zstore
    double precision, dimension(:), intent(in) :: xigll, yigll, zigll
    integer, dimension(:), intent(in) :: anchor_iax, anchor_iay, anchor_iaz
    logical,intent(in) :: POINT_CAN_BE_BURIED

    ! OUT
    integer,intent(out) :: ispec_selected
    double precision,intent(out) :: xi,eta,gamma
    double precision, intent(out) :: xix, xiy, xiz
    double precision, intent(out) :: etax, etay, etaz
    double precision, intent(out) :: gammax, gammay, gammaz
    double precision,intent(out) :: x,y,z
    double precision,intent(out) :: distmin_not_squared

    ! LOCAL parameters
    integer :: ix_initial_guess,iy_initial_guess,iz_initial_guess
    integer :: ispec,i,j,k
    integer(kind=8) :: iglob

    double precision :: lat,lon
    double precision :: distmin_squared,dist_squared

    logical :: target_located

    ! nodes search
    double precision,dimension(3) :: xyz_target
    double precision :: dist_min
    integer :: inode_min

    ! brute-force search for closest element
    ! if set to .false. (default), a kd-tree search for initial guess element is used
    logical,parameter :: DO_BRUTE_FORCE_SEARCH = .false.

    ! looks for closer estimates in neighbor elements if needed
    logical,parameter :: DO_ADJACENT_SEARCH = .true.

    ! logical, parameter :: DEBUG = .true.

    double precision, dimension(:,:), allocatable :: xyz_midpoints

    ! Get midpoints
    allocate(xyz_midpoints(nspec,3))
    xyz_midpoints(:,1) = xstore(ibool(MIDX, MIDY, MIDZ, :))
    xyz_midpoints(:,2) = ystore(ibool(MIDX, MIDY, MIDZ, :))
    xyz_midpoints(:,3) = zstore(ibool(MIDX, MIDY, MIDZ, :))

    if (DEBUG) write(IMAIN,*) "xyz midpoints"
    if (DEBUG) write(IMAIN,*) xyz_midpoints

    ! set distance to huge initial value
    distmin_squared = HUGEVAL

    ! initializes located target
    ! if we have not located a target element, the receiver is not in this slice
    ! therefore use first element only for fictitious iterative search
    ispec_selected = 1
    ix_initial_guess = MIDX
    iy_initial_guess = MIDY
    iz_initial_guess = MIDZ

    ! point in this slice

    ! kdtree search
    xyz_target(1) = x_target
    xyz_target(2) = y_target
    xyz_target(3) = z_target

    ! finds closest point (inside GLL points) in this chunk
    call kdtree_find_nearest_neighbor(xyz_target,inode_min,dist_min)

    ! debug
    !print *,'kd-tree found location :',inode_min

    ! sets element
    ispec_selected = inode_min

    ! loops over GLL points in this element to get (i,j,k) for initial guess
    do k = 2,NGLLZ-1
      do j = 2,NGLLY-1
        do i = 2,NGLLX-1
          iglob = ibool(i,j,k,ispec_selected)

          dist_squared = (x_target - dble(xstore(iglob)))*(x_target - dble(xstore(iglob))) &
                      + (y_target - dble(ystore(iglob)))*(y_target - dble(ystore(iglob))) &
                      + (z_target - dble(zstore(iglob)))*(z_target - dble(zstore(iglob)))

          if (DEBUG) write(IMAIN,*) "(i,j,k) dist2 (",i,j,k, ")", dist_squared
          ! take this point if it is closer to the receiver
          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          if (dist_squared < distmin_squared) then
            distmin_squared = dist_squared
            ix_initial_guess = i
            iy_initial_guess = j
            iz_initial_guess = k
          endif
        enddo
      enddo
    enddo

    if (DEBUG) write(IMAIN,*) "Ispec selected"
    if (DEBUG) write(IMAIN,*) ispec_selected
    if (DEBUG) write(IMAIN,*) "ijk first guess"
    if (DEBUG) write(IMAIN,*) ix_initial_guess, iy_initial_guess, iz_initial_guess

    ! ****************************************
    ! find the best (xi,eta,gamma)
    ! ****************************************


    !deprecated, but might be useful if one wishes to get an exact GLL point location:
    !
    !      ! for point sources, the location will be exactly at a GLL point
    !      ! otherwise this tries to find best location
    !
    !      if (USE_FORCE_POINT_SOURCE) then
    !        ! store xi,eta,gamma and x,y,z of point found
    !        ! note: they have range [1.0d0,NGLLX/Y/Z], used for point sources
    !        !          see e.g. in compute_add_sources.f90
    !        xi_subset(isource_in_this_subset) = dble(ix_initial_guess)
    !        eta_subset(isource_in_this_subset) = dble(iy_initial_guess)
    !        gamma_subset(isource_in_this_subset) = dble(iz_initial_guess)
    !
    !        iglob = ibool(ix_initial_guess,iy_initial_guess, &
    !            iz_initial_guess,ispec_selected)
    !        xyz_found_subset(1,isource_in_this_subset) = xstore(iglob)
    !        xyz_found_subset(2,isource_in_this_subset) = ystore(iglob)
    !        xyz_found_subset(3,isource_in_this_subset) = zstore(iglob)
    !
    !        ! compute final distance between asked and found (converted to km)
    !        final_distance_subset(isource_in_this_subset) = &
    !          dsqrt((x_target-xyz_found_subset(1,isource_in_this_subset))**2 + &
    !                (y_target-xyz_found_subset(2,isource_in_this_subset))**2 + &
    !                (z_target-xyz_found_subset(3,isource_in_this_subset))**2)*R_PLANET/1000.d0
    !      endif

    ! gets xi/eta/gamma and corresponding x/y/z coordinates
    call find_local_coordinates(x_target,y_target,z_target,xi,eta,gamma,x,y,z, &
                                xix, xiy, xiz, etax, etay, etaz, gammax, gammay, gammaz, &
                                ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess, &
                                nspec,NGLLX,NGLLY,NGLLZ, &
                                xstore,ystore,zstore, &
                                ibool, anchor_iax,anchor_iay,anchor_iaz, &
                                xigll,yigll,zigll, &
                                POINT_CAN_BE_BURIED)


    ! loops over neighbors and try to find better location
    if (DO_ADJACENT_SEARCH) then
      ! checks if position lies on an element boundary
      if (abs(xi) > 1.099d0 .or. abs(eta) > 1.099d0 .or. abs(gamma) > 1.099d0) then
        ! searches for better position in neighboring elements
        call find_best_neighbor(&
            x_target,y_target,z_target,xi,eta,gamma,x,y,z,ispec_selected,distmin_squared, &
            xix, xiy, xiz, etax, etay, etaz, gammax, gammay, gammaz, &
            nspec,NGLLX,NGLLY,NGLLZ, MIDX, MIDY, MIDZ,&
            xstore,ystore,zstore,xadj,adjncy, &
            ibool,anchor_iax,anchor_iay,anchor_iaz, &
            xigll,yigll,zigll, &
            POINT_CAN_BE_BURIED)
      endif
    endif ! DO_ADJACENT_SEARCH

    ! compute final distance between asked and found (converted to km)
    distmin_not_squared = dsqrt((x_target-x)*(x_target-x) + &
                                (y_target-y)*(y_target-y) + &
                                (z_target-z)*(z_target-z))*R_PLANET_KM

    !debug: sync
    !  endif ! iproc
    !  call synchronize_all()
    !  enddo

  end subroutine locate_point

!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_local_coordinates(x_target,y_target,z_target,xi,eta,gamma,x,y,z, &
                                    xix, xiy, xiz, etax, etay, etaz, gammax, gammay, gammaz, &
                                    ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess, &
                                    nspec,NGLLX,NGLLY,NGLLZ,xstore,ystore,zstore,ibool, &
                                    anchor_iax,anchor_iay,anchor_iaz,xigll,yigll,zigll, &
                                    POINT_CAN_BE_BURIED)

  use constants, only: NGNOD,HUGEVAL,NUM_ITER, IMAIN, DEBUG
  use jacobian, only: recompute_jacobian
  implicit none

  double precision,intent(in) :: x_target,y_target,z_target
  double precision,intent(out) :: xi,eta,gamma
  double precision,intent(out) :: x,y,z

  integer,intent(in) :: ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess

  integer, intent(in) :: nspec, NGLLX, NGLLY, NGLLZ
  integer, dimension(:), intent(in) :: anchor_iax, anchor_iay, anchor_iaz
  integer(kind=8), dimension(:,:,:,:), intent(in) :: ibool
  real(kind=4), dimension(:), intent(in) :: xstore, ystore, zstore
  double precision, dimension(:), intent(in) :: xigll, yigll, zigll


  logical,intent(in) :: POINT_CAN_BE_BURIED

  ! local parameters
  integer :: ia,iter_loop
  integer(kind=8) :: iglob

  ! coordinates of the control points of the surface element
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  double precision :: dx,dy,dz,dx_min,dy_min,dz_min,d_min_sq
  double precision :: dxi,deta,dgamma
  double precision, intent(inout) :: xix,xiy,xiz
  double precision, intent(inout) :: etax,etay,etaz
  double precision, intent(inout) :: gammax,gammay,gammaz

  ! define coordinates of the control points of the element
  do ia = 1,NGNOD
    iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec_selected)
    xelm(ia) = dble(xstore(iglob))
    yelm(ia) = dble(ystore(iglob))
    zelm(ia) = dble(zstore(iglob))
  enddo

  if (DEBUG) write (IMAIN,*) "Anchors"
  if (DEBUG) write (IMAIN,*) xelm
  if (DEBUG) write (IMAIN,*) yelm
  if (DEBUG) write (IMAIN,*) zelm

  ! use initial guess in xi and eta
  if (DEBUG) write (IMAIN,*) "xigll"
  if (DEBUG) write (IMAIN,*) xigll
  xi = xigll(ix_initial_guess)
  eta = yigll(iy_initial_guess)
  gamma = zigll(iz_initial_guess)

  if (DEBUG) write (IMAIN,*) 'initial guess'
  if (DEBUG) write (IMAIN,*) xi, eta, gamma
  if (DEBUG) write (IMAIN,*)

  ! impose receiver exactly at the surface
  if (.not. POINT_CAN_BE_BURIED) gamma = 1.d0

  d_min_sq = HUGEVAL
  dx_min = HUGEVAL
  dy_min = HUGEVAL
  dz_min = HUGEVAL

  ! iterate to solve the non linear system
  if (DEBUG) write (IMAIN,*) "START", x_target, y_target, z_target
  if (DEBUG) write (IMAIN,*) "     ", xi, eta, gamma
  if (DEBUG) write (IMAIN,*)

  if (DEBUG) write (IMAIN,*) "START   x ", x_target, y_target, z_target
  if (DEBUG) write (IMAIN,*) "       xi ", xi, eta, gamma
  if (DEBUG) write (IMAIN,*) "      xix ", xix,xiy,xiz
  if (DEBUG) write (IMAIN,*) "     etax ", etax,etay,etaz
  if (DEBUG) write (IMAIN,*) "   gammax ", gammax,gammay,gammaz

  ! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

    ! recompute Jacobian for the new point

    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                            xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

    ! compute distance to target location
    dx = - (x - x_target)
    dy = - (y - y_target)
    dz = - (z - z_target)

    if (DEBUG) write (IMAIN,*) '  iter ',iter_loop,'dx',sngl(dx),sngl(dx_min),'dy', &
                                sngl(dy),sngl(dy_min),'dz',sngl(dz),sngl(dz_min),d_min_sq

    ! compute increments
    if ((dx**2 + dy**2 + dz**2) < d_min_sq) then
      d_min_sq = dx**2 + dy**2 + dz**2
      dx_min = dx
      dy_min = dy
      dz_min = dz

      dxi = xix*dx + xiy*dy + xiz*dz
      deta = etax*dx + etay*dy + etaz*dz
      dgamma = gammax*dx + gammay*dy + gammaz*dz
    else
      ! new position is worse than old one, no change necessary
      dxi = 0.d0
      deta = 0.d0
      dgamma = 0.d0
    endif

    ! decreases step length if step is large
    if ((dxi**2 + deta**2 + dgamma**2) > 1.0d0) then
      dxi = dxi * 0.33333333333d0
      deta = deta * 0.33333333333d0
      dgamma = dgamma * 0.33333333333d0
    endif
    ! alternative: impose limit on increments (seems to result in slightly less accurate locations)
    !if (abs(dxi) > 0.3d0 ) dxi = sign(1.0d0,dxi)*0.3d0
    !if (abs(deta) > 0.3d0 ) deta = sign(1.0d0,deta)*0.3d0
    !if (abs(dgamma) > 0.3d0 ) dgamma = sign(1.0d0,dgamma)*0.3d0

    !debug
    !print *,'  dxi/..',(dxi**2 + deta**2 + dgamma**2),dxi,deta,dgamma

    ! update values
    xi = xi + dxi
    eta = eta + deta
    gamma = gamma + dgamma

    ! write (*,*) "    iter", iter_loop, x, y, z
    ! write (*,*) "                        ", xi, eta, gamma

    ! impose that we stay in that element
    ! (useful if user gives a receiver outside the mesh for instance)
    ! we can go slightly outside the [1,1] segment since with finite elements
    ! the polynomial solution is defined everywhere
    ! can be useful for convergence of iterative scheme with distorted elements
    if (xi > 1.10d0) xi = 1.10d0
    if (xi < -1.10d0) xi = -1.10d0
    if (eta > 1.10d0) eta = 1.10d0
    if (eta < -1.10d0) eta = -1.10d0
    if (gamma > 1.10d0) gamma = 1.10d0
    if (gamma < -1.10d0) gamma = -1.10d0

  ! end of non linear iterations
  enddo

  ! impose receiver exactly at the surface
  if (.not. POINT_CAN_BE_BURIED) gamma = 1.d0

  ! compute final coordinates of point found
  call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  end subroutine find_local_coordinates

!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_best_neighbor(&
    x_target,y_target,z_target,xi,eta,gamma,x,y,z,ispec_selected,distmin_squared, &
    xix, xiy, xiz, etax, etay, etaz, gammax, gammay, gammaz, &
    nspec,NGLLX,NGLLY,NGLLZ,MIDX,MIDY,MIDZ, &
    xstore,ystore,zstore, xadj, adjncy, &
    ibool, anchor_iax,anchor_iay,anchor_iaz, &
    xigll,yigll,zigll, &
    POINT_CAN_BE_BURIED)

  use constants, only: R_PLANET_KM
  use kdtree_search, only: kdtree_find_nearest_neighbor

  implicit none

  double precision,intent(in) :: x_target,y_target,z_target
  double precision,intent(inout) :: xi,eta,gamma
  double precision, intent(inout) :: xix,xiy,xiz
  double precision, intent(inout) :: etax,etay,etaz
  double precision, intent(inout) :: gammax,gammay,gammaz
  double precision,intent(inout) :: x,y,z

  integer,intent(inout) :: ispec_selected
  double precision,intent(inout) :: distmin_squared

  integer, intent(in) :: nspec, NGLLX, NGLLY, NGLLZ, MIDX, MIDY, MIDZ
  integer(kind=8), dimension(:,:,:,:), intent(in) :: ibool
  integer, dimension(:), intent(in) :: adjncy, xadj
  integer, dimension(:), intent(in) :: anchor_iax, anchor_iay, anchor_iaz
  real(kind=4), dimension(:), intent(in) :: xstore, ystore, zstore
  double precision, dimension(:), intent(in) :: xigll, yigll, zigll

  logical,intent(in) :: POINT_CAN_BE_BURIED

  ! local parameters
  integer :: ix_initial_guess,iy_initial_guess,iz_initial_guess
  integer :: ispec,i,j,k
  integer(kind=8) :: iglob

  double precision :: dist_squared
  double precision :: distmin_squared_guess

  ! nodes search
  double precision :: xi_n,eta_n,gamma_n,x_n,y_n,z_n ! neighbor position result

  ! neighbor elements
  integer,parameter :: MAX_NEIGHBORS = 50   ! maximum number of neighbors (around 37 should be sufficient for crust/mantle)
  integer :: index_neighbors(MAX_NEIGHBORS*MAX_NEIGHBORS) ! including neighbors of neighbors
  integer :: num_neighbors

  integer :: ii,jj,i_n,ientry,ispec_ref
  logical :: do_neighbor

  ! verbose output
  logical,parameter :: DEBUG = .true.

  ! best distance to target .. so far
  distmin_squared = (x_target - x)*(x_target - x) &
                  + (y_target - y)*(y_target - y) &
                  + (z_target - z)*(z_target - z)

  !debug
  if (DEBUG) print *,'neighbors: best guess ',ispec_selected,xi,eta,gamma,'distance',sngl(sqrt(distmin_squared)*R_PLANET_KM)

  ! fill neighbors arrays
  !
  ! note: we add direct neighbors plus neighbors of neighbors.
  !       for very coarse meshes, the initial location guesses especially around doubling layers can be poor such that we need
  !       to enlarge the search of neighboring elements.
  index_neighbors(:) = 0
  num_neighbors = 0
  do ii = 1,xadj(ispec_selected+1)-xadj(ispec_selected)
    ! get neighbor
    ientry = xadj(ispec_selected) + ii
    ispec_ref = adjncy(ientry)

    ! checks
    if (ispec_ref < 1 .or. ispec_ref > nspec) stop 'Invalid ispec index in locate point search'

    ! checks if exists already in list
    do_neighbor = .true.
    do i_n = 1,num_neighbors
      if (index_neighbors(i_n) == ispec_ref) then
        do_neighbor = .false.
        exit
      endif
    enddo

    ! adds to search elements
    if (do_neighbor) then
      num_neighbors = num_neighbors + 1
      index_neighbors(num_neighbors) = ispec_ref
    endif

    ! adds neighbors of neighbor
    do jj = 1,xadj(ispec_ref+1)-xadj(ispec_ref)
      ! get neighbor
      ientry = xadj(ispec_ref) + jj
      ispec = adjncy(ientry)

      ! checks
      if (ispec < 1 .or. ispec > nspec) stop 'Invalid ispec index in locate point search'

      ! checks if exists already in list
      do_neighbor = .true.
      do i_n = 1,num_neighbors
        if (index_neighbors(i_n) == ispec) then
          do_neighbor = .false.
          exit
        endif
      enddo

      ! adds to search elements
      if (do_neighbor) then
        num_neighbors = num_neighbors + 1
        index_neighbors(num_neighbors) = ispec
      endif
    enddo
  enddo

  ! loops over neighboring elements
  do i_n = 1,num_neighbors
    ispec = index_neighbors(i_n)

    ! note: the final position location can be still off if we start too far away.
    !       here we guess the best "inner" GLL point inside this search neighbor element
    !       to be the starting initial guess for finding the local coordinates.

    ! gets first guess as starting point
    ix_initial_guess = MIDX
    iy_initial_guess = MIDY
    iz_initial_guess = MIDZ
    iglob = ibool(MIDX,MIDY,MIDZ,ispec)

    distmin_squared_guess = (x_target - dble(xstore(iglob)))*(x_target - dble(xstore(iglob))) &
                          + (y_target - dble(ystore(iglob)))*(y_target - dble(ystore(iglob))) &
                          + (z_target - dble(zstore(iglob)))*(z_target - dble(zstore(iglob)))

    ! loop only on points inside the element
    ! exclude edges to ensure this point is not shared with other elements
    do k = 2,NGLLZ-1
      do j = 2,NGLLY-1
        do i = 2,NGLLX-1
          iglob = ibool(i,j,k,ispec)
          dist_squared = (x_target - dble(xstore(iglob)))*(x_target - dble(xstore(iglob))) &
                       + (y_target - dble(ystore(iglob)))*(y_target - dble(ystore(iglob))) &
                       + (z_target - dble(zstore(iglob)))*(z_target - dble(zstore(iglob)))
          !  keep this point if it is closer to the receiver
          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          if (dist_squared < distmin_squared_guess) then
            distmin_squared_guess = dist_squared
            ix_initial_guess = i
            iy_initial_guess = j
            iz_initial_guess = k
          endif
        enddo
      enddo
    enddo

    ! gets xi/eta/gamma and corresponding x/y/z coordinates
    call find_local_coordinates(x_target,y_target,z_target,xi_n,eta_n,gamma_n,x_n,y_n,z_n, &
                                xix, xiy, xiz, etax, etay, etaz, gammax, gammay, gammaz, &
                                ispec,ix_initial_guess,iy_initial_guess,iz_initial_guess, &
                                nspec,NGLLX,NGLLY,NGLLZ, &
                                xstore,ystore,zstore, &
                                ibool, &
                                anchor_iax, anchor_iay, anchor_iaz, &
                                xigll,yigll,zigll, &
                                POINT_CAN_BE_BURIED)

    ! final distance to target
    dist_squared = (x_target - x_n)*(x_target - x_n) &
                 + (y_target - y_n)*(y_target - y_n) &
                 + (z_target - z_n)*(z_target - z_n)

    ! debug
    if (DEBUG) print *,'  neighbor ',ispec,i_n,ientry,'ispec = ',ispec_selected,sngl(xi_n),sngl(eta_n),sngl(gamma_n), &
                       'distance',sngl(sqrt(dist_squared)*R_PLANET_KM),sngl(sqrt(distmin_squared)*R_PLANET_KM)

    ! takes this point if it is closer to the receiver
    ! (we compare squared distances instead of distances themselves to significantly speed up calculations)
    if (dist_squared < distmin_squared) then
      distmin_squared = dist_squared
      ! uses this as new location
      ispec_selected = ispec
      xi = xi_n
      eta = eta_n
      gamma = gamma_n
      x = x_n
      y = y_n
      z = z_n
    endif

    ! checks if position lies inside element (which usually means that located position is accurate)
    if (abs(xi) < 1.099d0 .and. abs(eta) < 1.099d0 .and. abs(gamma) < 1.099d0) exit

  enddo ! num_neighbors

  !debug
  if (DEBUG) print *,'neighbors: final ',ispec_selected,xi,eta,gamma,'distance',sqrt(distmin_squared)*R_PLANET_KM

  end subroutine find_best_neighbor

end module point_location