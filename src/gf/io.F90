submodule (gf) io
  implicit none
contains

  type(t_GF) module function read_GF(filename) result(GF)

    use hdf5
    use utils, only: throwerror
    use hdf5_utils, only: read_from_hdf5, get_dset_rank, get_dset_dims

    implicit none
    integer :: errorflag
    logical :: got
    character(len=100) :: errormessage

    ! file, group, attribute ids
    integer(hid_t)    :: file_id         ! file identifierd

    ! variable names
    character(len=65) :: filename ! input variable

    ! Random integers
    integer :: i,j,k

    ! ------ initialize hdf5 routines ----------------

    call h5open_f(errorflag)
    call throwerror(errorflag, " *** error initialising hdf routines")

    ! Opening file name
    call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, errorflag)
    call throwerror(errorflag, "Error when opening "//filename)

    ! Load Header variables
    call load_header(file_id, GF)

    ! Loading arrays
    call load_arrays(file_id, GF)

    ! ------ finalize routines -----------------------
    call h5fclose_f(file_id, errorflag)
    call throwerror(errorflag, "Error closing hdf5 routines")

    call h5close_f(errorflag)
    call throwerror(errorflag, "Error closing hdf5 routines")


  end function read_GF

  module subroutine print_GF(GF)

    type(t_GF), intent(in) :: GF

    ! Print format parameters
    CHARACTER(LEN=30) :: integerformat = "(A25 I20)"
    CHARACTER(LEN=30) :: realformat    = "(A25 F26.5)"
    CHARACTER(LEN=30) :: expformat     = "(A25 ES30.5)"

    write (*,*) "**************************************************************"
    write (*,*) "************************* HEADER *****************************"
    write (*,*) "**************************************************************"
    write (*,integerformat)   "do_adjacency_search:", GF%do_adjacency_search
    write (*,integerformat)   "nspec:",               GF%nspec
    write (*,integerformat)   "ngllx:",               GF%ngllx
    write (*,integerformat)   "nglly:",               GF%nglly
    write (*,integerformat)   "ngllz:",               GF%ngllz
    write (*,integerformat)   "midx:",                GF%midx
    write (*,integerformat)   "midy:",                GF%midy
    write (*,integerformat)   "midz:",                GF%midz
    write (*,integerformat)   "nsteps:",              GF%nsteps
    write (*,realformat)      "dt:",                  GF%dt
    write (*,realformat)      "hdur:",                GF%hdur
    write (*,realformat)      "tc:",                  GF%tc
    write (*,expformat)       "factor:",              GF%factor
    write (*,integerformat)   "topography:",          GF%topography

    if (GF%topography == 1) then
      write (*,realformat)      "resolution_topo_file:", GF%resolution_topo_file
      write (*,integerformat)   "nx_bathy:",             GF%nx_bathy
      write (*,integerformat)   "ny_bathy:",             GF%ny_bathy
    endif

    write (*,integerformat)   "ellipticity:",         GF%ellipticity

    write(*,*)
    write(*,*) "----------- Array info ----------------"
    write(*,*)
    write(*, integerformat)  "Number of elements:",   size(GF%ibool,4)
    write(*, integerformat)  "Number of GLL:",        size(GF%xyz,  1)
    write(*, integerformat)  "Ellipticity splines #", size(GF%rspl)
    write(*, integerformat)  "Number of stations",    size(GF%displacement, 1)
    write(*, integerformat)  "Number of components",  size(GF%displacement, 2)
    write(*,*)
    write(*,expformat)      "Displacment Min:", minval(GF%displacement)
    write(*,expformat)      "            Max:", maxval(GF%displacement)
    write(*,expformat)      "           Mean:", sum(GF%displacement)/size(GF%displacement)

    write (*,*) "**************************************************************"


  end subroutine print_GF


  module subroutine load_header(file_id, GF)

    use utils, only: throwerror
    use hdf5
    use hdf5_utils, only: read_from_hdf5, get_dset_dims

    ! file, group, attribute ids
    integer(hid_t), intent(in)    :: file_id

    ! GF structure to populate
    type(t_GF), intent(inout) :: GF

    ! Flags
    integer :: errorflag
    logical :: got

    ! Root name
    character(len=1)  :: name_root = "/"

    ! Network stations
    integer, dimension(1)              :: dims_stations, maxdims_stations
    integer, dimension(1)              :: dims_networks, maxdims_networks
    integer, dimension(1)              :: dims_burials, maxdims_burials
    integer, dimension(1)              :: dims_latitudes, maxdims_latitudes
    integer, dimension(1)              :: dims_longitudes, maxdims_longitudes

    character(len=30) :: name_do_adjacency_search
    character(len=30) :: name_ellipticity
    character(len=30) :: name_topography
    character(len=30) :: name_nsteps
    character(len=30) :: name_ngllx
    character(len=30) :: name_nglly
    character(len=30) :: name_ngllz
    character(len=30) :: name_nx_bathy
    character(len=30) :: name_ny_bathy
    character(len=30) :: name_factor
    character(len=30) :: name_dt
    character(len=30) :: name_hdur
    character(len=30) :: name_tc
    character(len=30) :: name_resolution_topo_file
    character(len=30) :: name_networks
    character(len=30) :: name_stations
    character(len=30) :: name_burials
    character(len=30) :: name_latitudes
    character(len=30) :: name_longitudes

    ! define variable names
    name_root                 = '/'
    name_do_adjacency_search  = name_root//'do_adjacency_search'
    name_ellipticity          = name_root//'ELLIPTICITY'
    name_topography           = name_root//'TOPOGRAPHY'
    name_nsteps               = name_root//'NSTEPS'
    name_ngllx                = name_root//'NGLLX'
    name_nglly                = name_root//'NGLLY'
    name_ngllz                = name_root//'NGLLZ'
    name_nx_bathy             = name_root//'NX_BATHY'
    name_ny_bathy             = name_root//'NY_BATHY'
    name_factor               = name_root//'FACTOR'
    name_dt                   = name_root//'DT'
    name_hdur                 = name_root//'HDUR'
    name_tc                   = name_root//'TC'
    name_resolution_topo_file = name_root//'RESOLUTION_TOPO_FILE'
    name_networks             = name_root//'Networks'
    name_stations             = name_root//'Stations'
    name_burials               = name_root//'burials'
    name_latitudes            = name_root//'latitudes'
    name_longitudes           = name_root//'longitudes'

    ! ------ Get scalar variables and flags
    call read_from_hdf5(GF%ellipticity, name_ellipticity, file_id, got, errorflag)
    call throwerror(errorflag, "Error reading "//name_ellipticity)


    call read_from_hdf5(GF%topography, name_topography, file_id, got, errorflag)
    call throwerror(errorflag, "Error reading "//name_topography)


    call read_from_hdf5(GF%do_adjacency_search, name_do_adjacency_search, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading "//name_do_adjacency_search)


    call read_from_hdf5(GF%ngllx, name_ngllx, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading "//name_ngllx)


    call read_from_hdf5(GF%nglly, name_nglly, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading "//name_nglly)


    call read_from_hdf5(GF%ngllz, name_ngllz, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading "//name_ngllz)



    call read_from_hdf5(GF%nsteps, name_nsteps, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading "//name_nsteps)

    call read_from_hdf5(GF%dt, name_dt, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading "//name_dt)


    call read_from_hdf5(GF%hdur, name_hdur, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading "//name_hdur)

    call read_from_hdf5(GF%tc, name_tc, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading "//name_tc)


    call read_from_hdf5(GF%factor, name_factor, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading "//name_factor)

    if (GF%topography == 1) then
      call read_from_hdf5(GF%nx_bathy, name_nx_bathy, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading "//name_nx_bathy)

      call read_from_hdf5(GF%ny_bathy, name_ny_bathy, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading "//name_ny_bathy)

      call read_from_hdf5(GF%resolution_topo_file, name_resolution_topo_file, &
                          file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading "//name_resolution_topo_file)
    endif

    ! Get dims of stations and networks
    call get_dset_dims(file_id, name_stations, dims_stations, maxdims_stations, got, errorflag)
    call throwerror(errorflag, "Error getting stations dims")
    call get_dset_dims(file_id, name_networks, dims_networks, maxdims_networks, got, errorflag)
    call throwerror(errorflag, "Error getting networks dims")
    call get_dset_dims(file_id, name_burials, dims_burials, maxdims_burials, got, errorflag)
    call throwerror(errorflag, "Error getting burials dims")
    call get_dset_dims(file_id, name_latitudes, dims_latitudes, maxdims_latitudes, got, errorflag)
    call throwerror(errorflag, "Error getting latitudes dims")
    call get_dset_dims(file_id, name_longitudes, dims_longitudes, maxdims_longitudes, got, errorflag)
    call throwerror(errorflag, "Error getting longitudes dims")

    ! Allocate stations and network names
    allocate(GF%networks(dims_networks(1)), GF%stations(dims_stations(1)))
    allocate(GF%burials(dims_burials(1)))
    allocate(GF%latitudes(dims_latitudes(1)), GF%longitudes(dims_longitudes(1)))

    ! Reading station and network names
    call read_from_hdf5(GF%stations, name_stations, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading stations")

    call read_from_hdf5(GF%networks, name_networks, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading networks")

    ! Reading latitude, longitude, and burial
    call read_from_hdf5(GF%latitudes, name_latitudes, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading latitudes")

    call read_from_hdf5(GF%longitudes, name_longitudes, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading longitudes")

    call read_from_hdf5(GF%burials, name_burials, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading burials")

    ! Finally define a midpoint number on the fly.
    GF%midx = GF%ngllx / 2 + 1
    GF%midy = GF%nglly / 2 + 1
    GF%midz = GF%ngllz / 2 + 1

  end subroutine load_header


  module subroutine load_arrays(file_id, GF)

    use constants, only: GAUSSALPHA, GAUSSBETA, IMAIN, DEBUG
    use gll_library, only: zwgljd
    use lagrange_poly, only: lagrange_any, lagrange_deriv_GLL
    use utils, only: throwerror
    use hdf5
    use hdf5_utils, only: read_from_hdf5, get_dset_dims

    implicit none

    ! file, group, attribute ids
    integer(hid_t), intent(in)    :: file_id

    ! GF structure to populate
    type(t_GF), intent(inout) :: GF

    ! Flags
    integer :: errorflag
    logical :: got

    ! iterators
    integer :: k1, k2, j1, j2, i1, i2

    ! Root name
    character(len=1)  :: name_root = "/"
    character(len=30) :: name_ellipticity_spline
    character(len=30) :: name_ellipticity_spline2
    character(len=30) :: name_rspl
    character(len=30) :: name_ibool
    character(len=30) :: name_bathy
    character(len=30) :: name_adjacency
    character(len=30) :: name_xadj
    character(len=30) :: name_displacement
    character(len=30) :: name_xyz

    ! Dimensions
    ! The actual dimensions are grabbed from file!
    ! So here the dimension only defines the rank of the arrays
    integer, dimension(1)              :: dims_ellipticity, maxdims_ellipticity
    integer, dimension(4)              :: dims_ibool, maxdims_ibool
    integer, dimension(2)              :: dims_xyz, maxdims_xyz
    integer, dimension(5)              :: dims_displacement, maxdims_displacement
    integer, dimension(2)              :: dims_bathy, maxdims_bathy
    integer, dimension(1)              :: dims_adjacency, maxdims_adjacency
    integer, dimension(1)              :: dims_xadj, maxdims_xadj

    ! Set names to be read
    name_rspl                 = name_root//'rspl'
    name_ellipticity_spline   = name_root//'ellipticity_spline'
    name_ellipticity_spline2  = name_root//'ellipticity_spline2'
    name_ibool                = name_root//'ibool'
    name_bathy                = name_root//'BATHY'
    name_adjacency            = name_root//'adjacency'
    name_xadj                 = name_root//'xadj'
    name_displacement         = name_root//'displacement'
    name_xyz                  = name_root//'xyz'

    ! integer, dimension(:,:,:,:)   :: temp_ibool
    ! integer, dimension(:,:,:,:,:) :: temp_diplacement
    ! integer, dimension(:,:,:,:)   :: temp_bathy

    ! ------ Get array dimensions ------------------

    ! Ibool
    call get_dset_dims(file_id, name_ibool, dims_ibool, maxdims_ibool, got, errorflag)
    call throwerror(errorflag, "Error getting ibool dims")


    ! Ellipticity Splines
    if (GF%ellipticity == 1) then
      call get_dset_dims(file_id, name_rspl, dims_ellipticity, maxdims_ellipticity, got, errorflag)
      call throwerror(errorflag, "Error getting ellipticity dims")

    endif

    ! Topography
    if (GF%topography == 1) then
      call get_dset_dims(file_id, name_bathy, dims_bathy, maxdims_bathy, got, errorflag)
      call throwerror(errorflag, "Error getting ellipticity dims")

    endif

    ! Topography
    call get_dset_dims(file_id, name_displacement, &
          dims_displacement, maxdims_displacement, got, errorflag)
    call throwerror(errorflag, "Error getting displacement dims")

    call get_dset_dims(file_id, name_xyz, &
          dims_xyz, maxdims_xyz, got, errorflag)
    call throwerror(errorflag, "Error getting xyz dims")

    call get_dset_dims(file_id, name_adjacency, &
          dims_adjacency, maxdims_adjacency, got, errorflag)
    call throwerror(errorflag, "Error getting adjacency dims")

    call get_dset_dims(file_id, name_xadj, &
          dims_xadj, maxdims_xadj, got, errorflag)
    call throwerror(errorflag, "Error getting xadj dims")


    ! ------ Allocate arrays ------------------------
    ! GLL interpolation values
    allocate(GF%xigll(GF%ngllx), GF%yigll(GF%nglly), GF%zigll(GF%ngllz), stat=errorflag)
    call throwerror(errorflag, "Error allocating interpolation values")

    ! GLL interpolation values
    allocate(GF%wxgll(GF%ngllx), GF%wygll(GF%nglly), GF%wzgll(GF%ngllz), stat=errorflag)
    call throwerror(errorflag, "Error allocating interpolation weights")

    ! GLL points and weights
    call zwgljd(GF%xigll(:),GF%wxgll(:),GF%ngllx,GAUSSALPHA,GAUSSBETA)
    call zwgljd(GF%yigll(:),GF%wygll(:),GF%nglly,GAUSSALPHA,GAUSSBETA)
    call zwgljd(GF%zigll(:),GF%wzgll(:),GF%ngllz,GAUSSALPHA,GAUSSBETA)

    if (GF%do_adjacency_search == 1) then
      allocate(GF%adjacency(dims_adjacency(1)), stat=errorflag)
      call throwerror(errorflag, "Error allocating adjacency vector")
      allocate(GF%xadj(dims_xadj(1)), stat=errorflag)
      call throwerror(errorflag, "Error allocating adjacency offsets")

      if (DEBUG) write(IMAIN,*) "dims adjacency: ", dims_adjacency(1)
      if (DEBUG) write(IMAIN,*) "dims xadj:      ", dims_xadj(1)
    endif

    if (GF%ellipticity == 1) then
      ! Spline radii
      allocate(GF%rspl(dims_ellipticity(1)), stat=errorflag)
      call throwerror(errorflag, "Error allocating ell. spline radii")

      ! Spline 1
      allocate(GF%ellipticity_spline(dims_ellipticity(1)), stat=errorflag)
      call throwerror(errorflag, "Error allocating ell. spline 1")

      ! Spline 2
      allocate(GF%ellipticity_spline2(dims_ellipticity(1)), stat=errorflag)
      call throwerror(errorflag, "Error allocating ell. spline 2")
    endif

    if (GF%topography == 1) then
      ! Topography
      allocate(GF%bathy(GF%nx_bathy, GF%ny_bathy), stat=errorflag)
      call throwerror(errorflag, "Error allocating topography")
    endif

    allocate(GF%ibool(dims_ibool(1), dims_ibool(2), dims_ibool(3), dims_ibool(4)), &
      stat=errorflag)
    call throwerror(errorflag, "Error allocating ibool")
    allocate(GF%displacement(dims_displacement(1), &
      dims_displacement(2), &
      dims_displacement(3), &
      dims_displacement(4), &
      dims_displacement(5)), &
      stat=errorflag)
    call throwerror(errorflag, "Error allocating displacement")

    allocate(GF%xyz(dims_xyz(1), dims_xyz(2)), stat=errorflag)
    call throwerror(errorflag, "Error allocating xyz coordinates")

    ! All is allocated let's read.

    call read_from_hdf5(GF%ibool, name_ibool, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading ibool")
    if (DEBUG) write(IMAIN,*) "ibool loaded."

    ! Add 1 to the indexing array
    GF%ibool(:,:,:,:) = GF%ibool(:,:,:,:) + 1

    call read_from_hdf5(GF%displacement, name_displacement, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading displacement")

    if (DEBUG) write(IMAIN,*) "post read min ", minval(GF%displacement)
    if (DEBUG) write(IMAIN,*) "post read max ", maxval(GF%displacement)

    call read_from_hdf5(GF%xyz, name_xyz, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading xyz")
    if (DEBUG) write(IMAIN,*) "XYZ loaded."

    if (GF%do_adjacency_search == 1) then

      call read_from_hdf5(GF%adjacency, name_adjacency, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading  adjacency")
      write(IMAIN,*) "adjacency loaded."

      call read_from_hdf5(GF%xadj, name_xadj, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading xadj")
      write(IMAIN,*) "xadj loaded."

    endif


    if (GF%ellipticity == 1) then
      call read_from_hdf5(GF%rspl, name_rspl, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading ellipticity spline")
      write(IMAIN,*) "rspl loaded."

      call read_from_hdf5(GF%ellipticity_spline, name_ellipticity_spline, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading ellipticity spline")
      write(IMAIN,*) "spline1 loaded."

      call read_from_hdf5(GF%ellipticity_spline2, name_ellipticity_spline2, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading ellipticity spline 2")
      write(IMAIN,*) "spline2 loaded."
    endif

    if (GF%topography == 1) then
      call read_from_hdf5(GF%bathy, name_bathy, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading bathymetry")
      write(IMAIN,*) "topobathy loaded."
    endif

    GF%nspec = dims_ibool(4)

  end subroutine load_arrays


  module subroutine free_GF(GF)

    ! GF structure to populate
    type(t_GF), intent(inout) :: GF

    if (allocated(GF%ibool)) then
      deallocate(GF%ibool)
    endif

    if (allocated(GF%rspl)) then
      deallocate(GF%rspl)
    endif

    if (allocated(GF%ellipticity_spline)) then
      deallocate(GF%ellipticity_spline)
    endif

    if (allocated(GF%ellipticity_spline2)) then
      deallocate(GF%ellipticity_spline2)
    endif

    if (allocated(GF%bathy)) then
      deallocate(GF%bathy)
    endif

    if (allocated(GF%displacement)) then
      deallocate(GF%displacement)
    endif

    if (allocated(GF%ibool)) then
      deallocate(GF%ibool)
    endif

    if (allocated(GF%adjacency)) then
      deallocate(GF%adjacency)
    endif

    if (allocated(GF%xadj)) then
      deallocate(GF%xadj)
    endif

    if (allocated(GF%xigll)) then
      deallocate(GF%xigll)
    endif

    if (allocated(GF%yigll)) then
      deallocate(GF%yigll)
    endif

    if (allocated(GF%zigll)) then
      deallocate(GF%zigll)
    endif

    if (allocated(GF%wxgll)) then
      deallocate(GF%wxgll)
    endif

    if (allocated(GF%wygll)) then
      deallocate(GF%wygll)
    endif

    if (allocated(GF%wzgll)) then
      deallocate(GF%wzgll)
    endif

  end subroutine free_GF

  !=====================================================================

end submodule io