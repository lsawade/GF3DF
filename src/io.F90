
module io

  private
  public :: read_GF, print_GF, free_GF, write_output_SAC

contains


  type(t_GF) function read_GF(filename) result(GF)

    use hdf5
    use utils, only: throwerror
    use ctypes, only: t_GF
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


  subroutine print_GF(GF)

    use ctypes, only: t_GF
    type(t_GF) :: GF

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

    write (*,*)
    write (*,*) "----------- Array info ----------------"
    write (*,*)
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


  subroutine load_header(file_id, GF)

    use ctypes, only: t_GF
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


    ! Finally define a midpoint number on the fly.
    GF%midx = GF%ngllx / 2 + 1
    GF%midy = GF%nglly / 2 + 1
    GF%midz = GF%ngllz / 2 + 1

  end subroutine load_header


  subroutine load_arrays(file_id, GF)

    use constants, only: GAUSSALPHA, GAUSSBETA
    use gll_library, only: zwgljd
    use lagrange_poly, only: lagrange_any, lagrange_deriv_GLL
    use ctypes, only: t_GF
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

    if (GF%topography == 1) then
      call get_dset_dims(file_id, name_bathy, dims_bathy, maxdims_bathy, got, errorflag)
      call throwerror(errorflag, "Error getting ellipticity dims")

    endif

    call get_dset_dims(file_id, name_displacement, &
          dims_displacement, maxdims_displacement, got, errorflag)
    call throwerror(errorflag, "Error getting displacement dims")


    call get_dset_dims(file_id, name_xyz, &
          dims_xyz, maxdims_xyz, got, errorflag)
    call throwerror(errorflag, "Error getting xyz dims")


    ! ------ Allocate arrays ------------------------


    ! GLL interpolation values
    allocate(GF%xigll(GF%ngllx), GF%yigll(GF%nglly), GF%zigll(GF%ngllz))
    allocate(GF%wxgll(GF%ngllx), GF%wygll(GF%nglly), GF%wzgll(GF%ngllz))

    ! GLL points and weights
    call zwgljd(GF%xigll(:),GF%wxgll(:),GF%ngllx,GAUSSALPHA,GAUSSBETA)
    call zwgljd(GF%yigll(:),GF%wygll(:),GF%nglly,GAUSSALPHA,GAUSSBETA)
    call zwgljd(GF%zigll(:),GF%wzgll(:),GF%ngllz,GAUSSALPHA,GAUSSBETA)

    if (GF%do_adjacency_search == 1) then
      allocate(GF%adjacency(dims_adjacency(1)))
      allocate(GF%xadj(dims_xadj(1)))
    endif

    if (GF%ellipticity == 1) then
      allocate(GF%rspl(dims_ellipticity(1)))
      allocate(GF%ellipticity_spline(dims_ellipticity(1)))
      allocate(GF%ellipticity_spline2(dims_ellipticity(1)))
    endif


    if (GF%topography == 1) then
      allocate(GF%bathy(GF%nx_bathy, GF%ny_bathy))
    endif
    write(*,*) "check point"

    allocate(GF%ibool(dims_ibool(1), dims_ibool(2), dims_ibool(3), dims_ibool(4)))
    allocate(GF%displacement(dims_displacement(1), &
      dims_displacement(2), &
      dims_displacement(3), &
      dims_displacement(4), &
      dims_displacement(5)))

    allocate(GF%xyz(dims_xyz(1), dims_xyz(2)))

    call read_from_hdf5(GF%ibool, name_ibool, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading ibool")

    ! Add 1 to the indexing array
    GF%ibool(:,:,:,:) = GF%ibool(:,:,:,:) + 1

    call read_from_hdf5(GF%displacement, name_displacement, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading displacement")

    write(*,*) "post read min ", minval(GF%displacement)
    write(*,*) "post read max ", maxval(GF%displacement)

    call read_from_hdf5(GF%xyz, name_xyz, file_id, got, errorflag)
    call throwerror(errorflag, "Error Reading xyz")

    if (GF%do_adjacency_search == 1) then

      call read_from_hdf5(GF%adjacency, name_adjacency, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading  adjacency")

      call read_from_hdf5(GF%xadj, name_xadj, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading xadj")

    endif


    if (GF%ellipticity == 1) then
      call read_from_hdf5(GF%rspl, name_rspl, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading ellipticity spline")

      call read_from_hdf5(GF%ellipticity_spline, name_ellipticity_spline, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading ellipticity spline")

      call read_from_hdf5(GF%ellipticity_spline2, name_ellipticity_spline2, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading ellipticity spline 2")
    endif

    if (GF%topography == 1) then
      call read_from_hdf5(GF%bathy, name_bathy, file_id, got, errorflag)
      call throwerror(errorflag, "Error Reading bathymetry")
    endif

    GF%nspec = dims_ibool(4)

  end subroutine load_arrays


  subroutine free_GF(GF)

    use ctypes, only: t_GF

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


  subroutine write_output_SAC(&
    seismograms, &
    iorientation, &
    sisname,&
    chn, &
    yr, &
    jda, &
    ho, &
    mi, &
    sec, &
    tshift_src, &
    t_shift, &
    event_name, &
    cmt_lat, &
    cmt_lon, &
    cmt_depth, &
    cmt_hdur, &
    station_name, &
    network_name, &
    stlat, &
    stlon, &
    stele, &
    stbur, &
    DT, &
    t0, &
    NT, &
    OUTPUT_SEISMOS_SAC_ALPHANUM, &
    OUTPUT_SEISMOS_SAC_BINARY, &
    MODEL, &
    OUTPUT_DIR)

    ! SAC headers have new format
    ! by Ebru

    use constants, only: IOUT_SAC, MAX_STRING_LEN
    use calendar, only: is_leap_year
    ! use binary_c_io, only: open_file_create, open_file_append, write_real, &
    !                        write_integer, write_character, write_n_real, &
    !                        close_file

    implicit none

    ! IN
    double precision, dimension(NT), intent(in) :: seismograms
    double precision, intent(in) :: t0, DT

    logical, intent(in) :: OUTPUT_SEISMOS_SAC_ALPHANUM
    logical, intent(in) :: OUTPUT_SEISMOS_SAC_BINARY
    character(len=MAX_STRING_LEN), intent(in) :: MODEL
    character(len=MAX_STRING_LEN), intent(in) :: OUTPUT_DIR

    ! Station parameters
    character(len=MAX_STRING_LEN), intent(in) :: sisname
    character(len=8), intent(in) :: network_name
    character(len=32), intent(in) :: station_name
    character(len=4), intent(in) :: chn
    integer, intent(in) :: iorientation
    double precision :: stlat, stlon, stele, stbur
    integer(kind=8) :: NT

    ! Event parameters
    character(len=*) :: event_name
    double precision :: cmt_lat , cmt_lon, cmt_depth, cmt_hdur
    double precision :: sec, t_shift, tshift_src
    integer :: yr, jda, ho, mi


    ! local parameters
    double precision :: btime
    integer :: seismo_offset
    real, dimension(NT) :: tmp
    integer :: time_sec,isample
    character(len=MAX_STRING_LEN) :: sisname_2

    ! SAC header variables
    real :: DELTA
    real :: DEPMIN
    real :: DEPMAX
    real :: SCALE_F
    real :: ODELTA
    real :: B,E,O,A
    real :: STLA,STLO,STEL,STDP
    real :: EVLA,EVLO,EVEL,EVDP
    real :: MAG,DIST,AZ,BAZ,GCARC
    real :: DEPMEN
    real :: USER0,USER1,USER2,USER3,USER4
    real :: CMPAZ,CMPINC

    integer :: NZYEAR,NZJDAY,NZHOUR,NZMIN,NZSEC
    integer :: NZMSEC,NVHDR,NORID,NEVID
    ! NUMBER of POINTS:
    integer :: NPTS
    integer :: IFTYPE,IMAGTYP
    integer :: IDEP
    integer :: IZTYPE
    integer :: IEVTYP
    integer :: IQUAL
    integer :: ISYNTH
    ! permission flags:
    integer :: LEVEN
    integer :: LPSPOL
    integer :: LOVROK
    integer :: LCALDA

    character(len=8) :: KSTNM
    character(len=16) :: KEVNM
    character(len=8) :: KCMPNM
    character(len=8) :: KNETWK
    character(len=8) :: KHOLE
    character(len=8) :: KUSER0,KUSER1,KUSER2
    character(len=8), parameter :: str_undef='-12345  '

    real :: UNUSED   ! header fields unused by SAC
    real :: undef    ! undefined values
    real :: INTERNAL ! SAC internal variables, always leave undefined
    real :: BYSAC
    ! end SAC header variables

    double precision :: shortest_period
    double precision :: value1,value2, value3,value4,value5

    integer :: imodulo_5


    !----------------------------------------------------------------

  !######################## SAC Alphanumeric Seismos ############################
  !
  ! written by Markus Treml and Bernhard Schuberth, Dept. for Earth and Environ-
  ! mental Sciences, Ludwig-Maximilians-University Munich, Germany
  !
  ! some words about SAC timing:
  !==============================
  !
  !NPTS,DELTA,B,E:
  ! These define the timing of the seismogram. E is calculated by sac. So, say
  ! you have 100 NPTS, a DELTA of 0.5, and set B to 0, E should be 50.
  ! Likewise setting B to -50 gives an E of 0.  Cutting basically cuts out points
  ! between the two times you designate based on these values.
  !KZTIME and KZDATE:
  ! Now things get funky.  KZTIME defines the exact time that the trace begins
  ! at. It has no affect on timing per se.  You'll really notice its effect if
  ! you read in two traces from different dates.

  ! Reference markers, (e.g. the o-marker) are not defined relative to this time,
  ! but rather to the begin time (B) of the seismo, so if you adjust B, you also
  ! need to adjust KZTIME to match. I would suggest experimenting with this until
  ! you understand it. It is a little non-intuitive until you see it for yourself.
  !
  !-----------------------------------------------------------------------------
  !
  ! This file is essentially the alphanumeric equivalent of the SAC binary data
  ! file. The header section is stored on the first 30 cards. This is followed
  ! by one or two data sections. The data is in 5G15.7 format.
  !----------------------------------------------------------------------
  !
  ! SAC header file format: https://ds.iris.edu/files/sac-manual/manual/file_format.html

    ! define certain default values

    ! unused or undefined values are set to '-12345.00'
    UNUSED   = -12345.00 ! header fields unused by SAC
    undef    = -12345.00 ! undefined values
    INTERNAL = -12345.00 ! SAC internal variables, always left undefined
    BYSAC    = -12345.00 ! values calculated by SAC from other variables
    !
    DELTA  = sngl(DT)    ! [REQUIRED]
    DEPMIN = BYSAC
    DEPMAX = BYSAC
    DEPMEN = BYSAC
    SCALE_F= 1000000000  ! factor for y-value, set to 10e9, so that values are in nm
    ODELTA = undef       ! increment from delta

    ! begin time
    btime = (seismo_offset)*DT - t0 + tshift_src

    B      = sngl(btime) ! [REQUIRED]
    E      = BYSAC       ! [REQUIRED]
    O      = 0  !
    A      = undef  !###

    !station values:
    STLA = sngl(stlat)
    STLO = sngl(stlon)
    STEL = sngl(stele)
    STDP = sngl(stbur)

    !event values (hypocenter):
    ! note: this writes out the CMT location, which might be different
    ! to the event location given in the first, PDE line
    EVLA   = sngl(cmt_lat)
    EVLO   = sngl(cmt_lon)
    EVEL   = undef  !not defined
    EVDP   = sngl(cmt_depth)

    ! SAC headers will have new format
    USER0  = sngl(cmt_hdur) !half duration from CMT file if not changed to t0 = 0.d0 (point source)

    ! USER1 and USER2 slots are used for the shortest and longest periods at which
    ! simulations are accurate, respectively.

    ! minimum period estimation
    shortest_period = 1/DT

    USER1  = sngl(shortest_period)
    USER2  = 500.0d0
    USER3  = undef
    USER4  = undef
    ! we remove any PDE information, since the simulation could also start
    ! with a "pure" CMT solution, without having any PDE info
    !
    !USER1  = sngl(t_shift) !time shift between PDE and CMT solutions
    !PDE location values (different from CMT location, usually):
    !USER2  = sngl(depth) !PDE depth
    !USER3  = sngl(elat) !PDE event latitude
    !USER4  = sngl(elon) !PDE event longitude
    !
    !cmt location values (different from hypocenter location, usually):
    ! USER0  = sngl(cmt_lat)
    ! USER1  = sngl(cmt_lon)
    !USER0  = sngl(elat)
    !USER1  = sngl(elon)
    !USER2  = sngl(depth)
    !USER3  = sngl(cmt_hdur) !half duration from CMT if not changed to t0 = 0.d0 (point source)

    ! Initialize values
    value1 = 0.d0
    value2 = 0.d0
    value3 = 0.d0
    value4 = 0.d0
    value5 = 0.d0

    ! it is not clear, which magnitude to write out:
    ! should it be
    !   body-wave-magnitude (Mb), surface-wave-magnitude (Ms), moment magnitude (Mw)
    !   or leave magnitude and use scalar moment (M0, but calculated by which convention, Harvard?)
    !
    ! it's confusing, and as a result, we will omit it.
    MAG     = undef
    IMAGTYP = int(undef)

    !MAG    = mb    !
    !IMAGTYP= 52    ! 52 = Mb? 55 = Mw!

    DIST   = BYSAC ! cause
    AZ     = BYSAC ! LCALDA
    BAZ    = BYSAC ! is
    GCARC  = BYSAC ! TRUE

    ! instrument orientation
    if (iorientation == 1) then !N
      CMPAZ  = 0.00
      CMPINC =90.00
    else if (iorientation == 2) then !E
      CMPAZ  =90.00
      CMPINC =90.00
    else if (iorientation == 3) then !Z
      CMPAZ  = 0.00
      CMPINC = 0.00
    endif
    !----------------end format G15.7--------

    ! date and time:
    NZYEAR = yr
    NZJDAY = jda
    NZHOUR = ho
    NZMIN  = mi

    ! adds time-shift to get the CMT time in the headers as origin time of events
    NZSEC  = int(sec+t_shift)
    NZMSEC = int((sec+t_shift-int(sec+t_shift))*1000)

    !NZSEC  =int(sec)
    !NZMSEC =int((sec-int(sec))*1000)

    ! Adjust event time and date after t_shift is added
    if (NZSEC >= 60) then
    time_sec = jda*24*3600 + ho*3600 + mi*60 + int(sec+t_shift)
    NZJDAY   = int(time_sec/(24*3600))
    NZHOUR   = int(mod(time_sec,24*3600)/3600)
    NZMIN    = int(mod(time_sec,3600)/60)
    NZSEC    = mod(time_sec,60)
    if (NZJDAY > 365 .and. .not. is_leap_year(NZYEAR)) then
        NZJDAY = mod(NZJDAY,365)
        NZYEAR = yr + 1
    else if (NZJDAY > 366 .and. is_leap_year(NZYEAR)) then
        NZJDAY = mod(NZJDAY,366)
        NZYEAR = yr + 1
    else if (NZJDAY == 366 .and. is_leap_year(NZYEAR)) then
        NZJDAY = 366
    endif
    endif


    NVHDR=6 ! SAC header version number. Current is 6

    ! CSS3.0 variables:
    NORID = int(undef) !origin ID
    NEVID = int(undef) !event  ID
    !NWVID =undef !waveform ID

    ! NUMBER of POINTS:
    NPTS = NT ! [REQUIRED]
    ! event type
    IFTYPE = 1 ! 1=ITIME, i.e. seismogram  [REQUIRED] # numbering system is
    IDEP   = 6 ! 6: displ/nm                          # quite strange, best

    IZTYPE = 11 !=origint reference time equivalent ! # by chnhdr and write
    IEVTYP = 40 !event type, 40: Earthquake           # alpha and check
    IQUAL  = int(undef) ! quality
    ISYNTH = int(undef) ! 1 real data, 2...n synth. flag
    ! permission flags:
    LEVEN  = 1 ! evenly spaced data [REQUIRED]
    LPSPOL = 1 ! ? pos. polarity of components (has to be TRUE for LCALDA=1)
    LOVROK = 1 ! 1: OK to overwrite file on disk
    LCALDA = 1 ! 1: calculate DIST, AZ, BAZ, and GCARC, 0: do nothing
    ! ------------------end format 5I10---------
    !
    !----------------------------------
    KSTNM  = station_name(1:8) ! A8

    ! writes out event id as event name
    KEVNM  = event_name(1:len_trim(event_name)) ! A16

    !if (NSOURCES == 1) then
    !  KEVNM  = ename(1:len_trim(ename))//'_syn'! A16
    !else
    !  KEVNM  = ename(1:len_trim(ename))//'_sFS'! A16
    !endif

    KCMPNM = chn(1:3)           ! 3A8
    KNETWK = network_name(1:6)  !  A6

    ! KHOLE slot represents SEED location IDs.
    ! Based on the IRIS convention, S1 and S3 are assigned to 1D and 3D seismograms, respectively.
    ! If a model is a combination of 1D and 3D models (e.g., 3D mantle with 1D crust), it will be considered as 3D.
    ! Currently, the decision is made based on model names given in Par_file assuming that
    ! all 1D model names start with "1D".
    ! Ebru, December 1, 2011

    KHOLE = 'S3'
    if (trim(MODEL(1:2)) == "1D") KHOLE = 'S1'

    ! indicates SEM synthetics
    KUSER0 = 'SY'          ! Network code assigned by IRIS for synthetic seismograms
    KUSER1 = 'SEM8.0.0'    ! code version 8.0
    KUSER2 = 'Tiger'       ! year of the tiger: Feb 01 2022 - Jan 21 2023
                          ! (chinese zodiac http://en.wikipedia.org/wiki/Chinese_zodiac :)

    !KUSER0 = 'PDE_LAT_'          !  A8
    !KUSER1 = 'PDE_LON_'          !  A8
    !KUSER2 = 'PDEDEPTH'          !  A8
    !----------------------------------

    if (OUTPUT_SEISMOS_SAC_ALPHANUM) then

      ! add .sacan (sac alphanumeric) extension to seismogram file name for SAC seismograms
      write(sisname_2,"('/',a,'.sacan')") trim(sisname)

      open(unit=IOUT_SAC,file=trim(OUTPUT_DIR)//trim(sisname_2), &
        status='unknown',action='write')

  ! Formats of alphanumerical SAC header fields
  510 format(5G15.7,5G15.7,5G15.7,5G15.7,5G15.7)
  520 format(5I10,5I10,5I10,5I10,5I10)
  530 format(A8,A16)
  540 format(A8,A8,A8)


      if (seismo_offset == 0) then
        !
        ! now write actual header:
        ! ------------------------
        !
        ! real variables:
        !                                 DELTA     DEPMIN   DEPMAX   SCALE   ODELTA
        !                                 B         E        O        A       INTERNAL
        !                                 T0        T1       T2       T3      T4
        !                                 T5        T6       T7       T8      T9
        !                                 F         RESP0    RESP1    RESP2   RESP3
        !                                 RESP4     RESP5    RESP6    RESP7   RESP8
        !                                 RESP9     STLA     STLO     STEL    STDP
        !                                 EVLA      EVLO     EVEL     EVDP    MAG
        !                                 USER0     USER1    USER2    USER3   USER4
        !                                 USER5     USER6    USER7    USER8   USER9
        !                                 DIST      AZ       BAZ      GCARC   INTERNAL
        !                                 INTERNAL  DEPMEN   CMPAZ    CMPINC  XMINIMUM
        !                                 XMAXIMUM  YMINIMUM YMAXIMUM ADJTM   UNUSED
        !
        write(IOUT_SAC,510) DELTA,    DEPMIN,  DEPMAX,  SCALE_F,  ODELTA
        write(IOUT_SAC,510) B,        E,       O,       A,      INTERNAL
        write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
        write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
        write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
        write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
        write(IOUT_SAC,510) undef,    STLA,    STLO,    STEL,   STDP
        write(IOUT_SAC,510) EVLA,     EVLO,    EVEL,    EVDP,   MAG
        write(IOUT_SAC,510) USER0,    USER1,   USER2,   USER3,  USER4
        write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
        write(IOUT_SAC,510) DIST,     AZ,      BAZ,     GCARC,  INTERNAL
        write(IOUT_SAC,510) INTERNAL, DEPMEN,  CMPAZ,   CMPINC, undef
        write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
        write(IOUT_SAC,510) UNUSED,   UNUSED,  UNUSED,  UNUSED, UNUSED
        !
        ! integer variables:
        !                                 NSPTS, NWFID, NXSIZE, NYSIZE, UNUSED
        !                                                                    IINST
        !                                 ISTREG IEVREG IEVTYP IQUAL ISYNTH
        !                                 IMAGTYP, IMAGSRC, UNUSED, UNUSED, UNUSED
        !
        write(IOUT_SAC,520) NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC
        write(IOUT_SAC,520) NZMSEC, NVHDR, NORID, NEVID, NPTS
        write(IOUT_SAC,520) int(undef),int(undef),int(undef),int(undef),int(undef)
        write(IOUT_SAC,520) IFTYPE, IDEP, IZTYPE, int(UNUSED), int(undef)
        write(IOUT_SAC,520) int(undef),int(undef),IEVTYP, int(undef), ISYNTH
        write(IOUT_SAC,520) IMAGTYP,int(undef),int(undef),int(undef),int(undef)
        write(IOUT_SAC,520) int(UNUSED), int(UNUSED), int(UNUSED), int(UNUSED), int(UNUSED)
        write(IOUT_SAC,520) LEVEN, LPSPOL, LOVROK, LCALDA, int(UNUSED)
        write(IOUT_SAC,530) KSTNM, KEVNM
        !
        ! character variables:
        !
        !                                   KHOLE    KO       KA
        !                                   KT0      KT1      KT2
        !                                   KT3      KT4      KT5
        !                                   KT6      KT7      KT8
        !                                   KT9      KF       KUSER0
        !                                   KUSER1     KUSER2       KCMPNM
        !                                   KNETWK   KDATRD   KINST
        !
        write(IOUT_SAC,540) KHOLE,'-12345  ','-12345  '
        write(IOUT_SAC,540) '-12345  ','-12345  ','-12345  '
        write(IOUT_SAC,540) '-12345  ','-12345  ','-12345  '
        write(IOUT_SAC,540) '-12345  ','-12345  ','-12345  '
        write(IOUT_SAC,540) '-12345  ','-12345  ',KUSER0
        write(IOUT_SAC,540)   KUSER1, KUSER2, KCMPNM
        write(IOUT_SAC,540)   KNETWK,'-12345  ','-12345  '
      endif

      ! now write data - with five values per row:
      ! ---------------
      imodulo_5 = mod(NT,5)
      if (imodulo_5 == 0) then
        ! five values per row
        do isample = 1+5,NT+1,5
          value1 = dble(seismograms(isample-5))
          value2 = dble(seismograms(isample-4))
          value3 = dble(seismograms(isample-3))
          value4 = dble(seismograms(isample-2))
          value5 = dble(seismograms(isample-1))
          write(IOUT_SAC,510) sngl(value1),sngl(value2),sngl(value3),sngl(value4),sngl(value5)
        enddo
      else
        ! five values per row as long as possible
        do isample = 1+5,(NT-imodulo_5)+1,5
          value1 = dble(seismograms(isample-5))
          value2 = dble(seismograms(isample-4))
          value3 = dble(seismograms(isample-3))
          value4 = dble(seismograms(isample-2))
          value5 = dble(seismograms(isample-1))
          write(IOUT_SAC,510) sngl(value1),sngl(value2),sngl(value3),sngl(value4),sngl(value5)
        enddo
        ! loads remaining values
        if (imodulo_5 >= 1) value1 = dble(seismograms(NT-imodulo_5+1))
        if (imodulo_5 >= 2) value2 = dble(seismograms(NT-imodulo_5+2))
        if (imodulo_5 >= 3) value3 = dble(seismograms(NT-imodulo_5+3))
        if (imodulo_5 >= 4) value4 = dble(seismograms(NT-imodulo_5+4))
        ! writes out last data line
        select case(imodulo_5)
        case (1)
          write(IOUT_SAC,510) sngl(value1)
        case (2)
          write(IOUT_SAC,510) sngl(value1),sngl(value2)
        case (3)
          write(IOUT_SAC,510) sngl(value1),sngl(value2),sngl(value3)
        case (4)
          write(IOUT_SAC,510) sngl(value1),sngl(value2),sngl(value3),sngl(value4)
        case default
          stop 'Error invalid SAC alphanumeric output'
        end select
      endif

    endif ! OUTPUT_SEISMOS_SAC_ALPHANUM

    ! For explanation on values set, see above (SAC ASCII)
    if (OUTPUT_SEISMOS_SAC_BINARY) then

      ! add .sac (sac binary) extension to seismogram file name for SAC seismograms
      write(sisname_2,"('/',a,'.sac')") trim(sisname)

      ! open binary file
      if (seismo_offset == 0) then
        call open_file_create(trim(OUTPUT_DIR)//trim(sisname_2)//char(0))
      else
        call open_file_append(trim(OUTPUT_DIR)//trim(sisname_2)//char(0))
      endif

      if (seismo_offset == 0) then
        ! write header variables

        ! write single precision header variables 1:70
        call write_real(DELTA)         !(1)
        call write_real(DEPMIN)        !(2)
        call write_real(DEPMAX)        !(3)
        call write_real(SCALE_F)       !(4)
        call write_real(ODELTA)        !(5)
        call write_real(B)             !(6)
        call write_real(E)             !(7)
        call write_real(O)             !(8)
        call write_real(A)             !(9)
        call write_real(INTERNAL)      !(10)
        call write_real(undef)          !(11)T0
        call write_real(undef)          !(12)T1
        call write_real(undef)          !(13)T2
        call write_real(undef)          !(14)T3
        call write_real(undef)          !(15)T4
        call write_real(undef)          !(16)T5
        call write_real(undef)          !(17)T6
        call write_real(undef)          !(18)T7
        call write_real(undef)          !(19)T8
        call write_real(undef)          !(20)T9
        call write_real(undef)          !(21)F
        call write_real(undef)          !(22)RESP0
        call write_real(undef)          !(23)RESP1
        call write_real(undef)          !(24)RESP2
        call write_real(undef)          !(25)RESP3
        call write_real(undef)          !(26)RESP4
        call write_real(undef)          !(27)RESP5
        call write_real(undef)          !(28)RESP6
        call write_real(undef)          !(29)RESP7
        call write_real(undef)          !(30)RESP8
        call write_real(undef)          !(31)RESP9
        call write_real(STLA)          !(32)
        call write_real(STLO)          !(33)
        call write_real(STEL)          !(34)
        call write_real(STDP)          !(35)
        call write_real(EVLA)          !(36)
        call write_real(EVLO)          !(37)
        call write_real(EVEL)          !(38)
        call write_real(EVDP)          !(39)
        call write_real(MAG)           !(40)
        call write_real(USER0)         !(41)USER0
        call write_real(USER1)         !(42)USER1
        call write_real(USER2)         !(43)USER2
        call write_real(undef)         !(44)USER3
        call write_real(undef)          !(45)USER4
        call write_real(undef)          !(46)USER5
        call write_real(undef)          !(47)USER6
        call write_real(undef)          !(48)USER7
        call write_real(undef)          !(49)USER8
        call write_real(undef)          !(50)USER9
        call write_real(DIST)          !(51)
        call write_real(AZ)            !(52)
        call write_real(BAZ)           !(53)
        call write_real(GCARC)         !(54)
        call write_real(INTERNAL)      !(55)
        call write_real(INTERNAL)      !(56)
        call write_real(DEPMEN)        !(57)
        call write_real(CMPAZ)         !(58)
        call write_real(CMPINC)        !(59)
        call write_real(undef)          !(60)XMINIMUM
        call write_real(undef)          !(61)XMAXIMUM
        call write_real(undef)          !(62)YMINIMUM
        call write_real(undef)          !(63)YMAXIMUM
        call write_real(undef)          !(64)
        call write_real(undef)          !(65)
        call write_real(undef)          !(66)
        call write_real(undef)          !(67)
        call write_real(undef)          !(68)
        call write_real(undef)          !(69)
        call write_real(undef)          !(70)

        ! write integer header variables 71:105
        call write_integer(NZYEAR)        !(71)
        call write_integer(NZJDAY)        !(72)
        call write_integer(NZHOUR)        !(73)
        call write_integer(NZMIN)         !(74)
        call write_integer(NZSEC)         !(75)
        call write_integer(NZMSEC)        !(76)
        call write_integer(NVHDR)         !(77)
        call write_integer(NORID)         !(78)
        call write_integer(NEVID)         !(79)
        call write_integer(NPTS)          !(80)
        call write_integer(int(undef))     !(81)UNUSED
        call write_integer(int(undef))     !(82)NWFID
        call write_integer(int(undef))     !(83)NXSIZE
        call write_integer(int(undef))     !(84)NYSIZE
        call write_integer(int(undef))     !(85)UNUSED
        call write_integer(IFTYPE)        !(86)
        call write_integer(IDEP)          !(87)
        call write_integer(IZTYPE)        !(88)
        call write_integer(int(undef))     !(89)UNUSED
        call write_integer(int(undef))     !(90)IINST
        call write_integer(int(undef))     !(91)ISTREG
        call write_integer(int(undef))     !(92)IEVREG
        call write_integer(IEVTYP)        !(93)
        call write_integer(IQUAL)         !(94)
        call write_integer(ISYNTH)        !(95)
        call write_integer(IMAGTYP)       !(96)
        call write_integer(int(undef))     !(97)IMAGSRC
        call write_integer(int(UNUSED))   !(98)
        call write_integer(int(UNUSED))   !(99)
        call write_integer(int(UNUSED))   !(100)
        call write_integer(int(UNUSED))   !(101)
        call write_integer(int(UNUSED))   !(102)
        call write_integer(int(UNUSED))   !(103)
        call write_integer(int(UNUSED))   !(104)
        call write_integer(int(UNUSED))   !(105)

        ! write logical header variables 106:110
        call write_integer(LEVEN)         !(106)
        call write_integer(LPSPOL)        !(107)
        call write_integer(LOVROK)        !(108)
        call write_integer(LCALDA)        !(109)
        call write_integer(int(UNUSED))   !(110)


        ! write character header variables 111:302
        call write_character(KSTNM,8)         !(111:118)
        call write_character(KEVNM,16)         !(119:134)
        call write_character(KHOLE,8)          !(135:142)KHOLE
        call write_character(str_undef,8)      !(143:150)KO
        call write_character(str_undef,8)      !(151:158)KA
        call write_character(str_undef,8)      !(159:166)KT0
        call write_character(str_undef,8)      !(167:174)KT1
        call write_character(str_undef,8)      !(175:182)KT2
        call write_character(str_undef,8)      !(183:190)KT3
        call write_character(str_undef,8)      !(191:198)KT4
        call write_character(str_undef,8)      !(199:206)KT5
        call write_character(str_undef,8)      !(207:214)KT6
        call write_character(str_undef,8)      !(215:222)KT7
        call write_character(str_undef,8)      !(223:230)KT8
        call write_character(str_undef,8)      !(231:238)KT9
        call write_character(str_undef,8)      !(239:246)KF
        call write_character(KUSER0,8)        !(247:254)
        call write_character(KUSER1,8)        !(255:262)
        call write_character(KUSER2,8)        !(263:270)
        call write_character(KCMPNM,8)        !(271:278)
        call write_character(KNETWK,8)        !(279:286)
        call write_character(str_undef,8)      !(287:294)KDATRD
        call write_character(str_undef,8)      !(295:302)KINST

      endif

      ! now write SAC time series to file
      ! BS BS write whole time series at once (hope to increase I/O performance
      ! compared to using a loop on it)
      tmp(1:NT) = real(seismograms(1:NT))
      call write_n_real(tmp(1:NT),NT)

      call close_file()

    endif ! OUTPUT_SEISMOS_SAC_BINARY

  end subroutine write_output_SAC


end module io