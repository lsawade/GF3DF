
program open_file

  use hdf5
!   include 'hdf5_fortran'
  use hdf5_utils

  implicit none
  integer :: errorflag
  logical :: got
  character(len=100) :: errormessage

  ! file, group, attribute ids
  integer(hid_t)    :: file_id         ! file identifier
  integer(hid_t)    :: root_id         ! root id
  integer(hid_t)    :: ellipticity_id  ! variables
  integer(hid_t)    :: factor_id       !    ...

  ! variable names
  character(len=65) :: filename ! input variable
  character(len=1)  :: root_name
  character(len=13) :: ellipticity_name
  character(len=8)  :: factor_name
  character(len=6)  :: ibool_name
  character(len=6)  :: ngllx_name
  character(len=6)  :: nglly_name
  character(len=6)  :: ngllz_name
  character(len=6)  :: nglob_name


  ! variable types
  integer(hsize_t), dimension(1) :: dim1
  integer,          dimension(1) :: ellipticity
  integer                        :: ngllx = 0
  integer                        :: nglly = 0
  integer                        :: ngllz = 0
  integer                        :: nglob = 0
  integer :: ellipticity_util
  double precision, dimension(1) :: factor
  integer, allocatable, dimension(:,:,:,:) :: ibool

  ! define variable names
  root_name = '/'
  ellipticity_name = root_name//'ELLIPTICITY'
  factor_name = root_name//'FACTOR'
  ibool_name = root_name//'ibool'
  ngllx_name = root_name//'NGLLX'
  nglly_name = root_name//'NGLLY'
  ngllz_name = root_name//'NGLLZ'
  nglob_name = root_name//'NGLOB'

  ! dimensions
  dim1(1) = 1

  ! ------ initialize hdf5 routines ----------------

  call h5open_f(errorflag)

  if (errorflag/=0) then
     errormessage=" *** error initialising hdf routines"
     return
  endif

  ! ------ open file hdf5 routines ----------------

  filename = '../single_element.h5' ! refer to line 5 of hdf5.js
  ! which gives the filename

  call h5fopen_f(filename, h5f_acc_rdonly_f, file_id, errorflag)

  if (errorflag/=0) then
     errormessage=" *** error opening hdf file"
     return
  endif

  ! ------ open root group -------------------------

  call h5gopen_f(file_id, root_name, root_id, errorflag)

  if (errorflag/=0) then
     errormessage=" *** error opening hdf file"
     return
  endif

  ! ------ read attritbute ------------------------

  call h5dopen_f(file_id, ellipticity_name, ellipticity_id, errorflag)

  if (errorflag/=0) then
     errormessage=" *** error opening total attribute "
     return
  endif

  call h5dread_f(ellipticity_id, h5t_native_integer, ellipticity, dim1, errorflag)

  if (errorflag/=0) then
     errormessage=" *** error opening total attribute "
     return
  endif

  write (*,*) "this is ellipticity: ", ellipticity

  call read_from_hdf5(ellipticity_util, ellipticity_name, ellipticity_id, got, errorflag)

  if (errorflag/=0) then
     errormessage=" *** error opening total attribute "
     return
  endif

  write (*,*) "this is ellipticity util: ", ellipticity_util

  call read_from_hdf5(ngllx, ngllx_name, file_id, got, errorflag)

  if (errorflag/=0) then
     write (*,*) " *** error opening ngllx"
  endif
  write (*,*) "got, ", got

  call read_from_hdf5(nglly, nglly_name, file_id, got, errorflag)

  if (errorflag/=0) then
     write (*,*) " *** error opening nglly"
  endif
    write (*,*) "got, ", got

  call read_from_hdf5(ngllz, ngllz_name, file_id, got, errorflag)

  if (errorflag/=0) then
     write (*,*) " *** error opening ngllz"
  endif

  write (*,*) "got, ", got
  write (*,*) "GLL structure (NGLLX, NGLLY, NGLLZ): (", ngllx, ",", nglly, ",", nglly,")"

!   call read_from_hdf5(nglob, nglob_name, file_id, got, errorflag)

!   if (errorflag/=0) then
!      write (*,*) " *** error opening nglob"
!   endif

  write (*,*) "got, ", got
  write (*,*) "NGLOB", nglob

!   allocate (ibool(ngllx, nglly, ngllz, nglob))
!   call read_from_hdf5(ibool, )

  ! ------ close file hdf5 -------------------------

  call h5fclose_f(file_id, errorflag)

  if (errorflag/=0) then
     errormessage=" *** error opening hdf file"
     return
  endif

  ! ------ finalize routines -----------------------

  call h5close_f(errorflag)

  if (errorflag /= 0) then
     errormessage=" *** error opening hdf file"
     return
  endif


 end program open_file
