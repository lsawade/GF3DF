module hdf5_utils

  use hdf5
  use iso_c_binding, only:c_loc
  implicit none

  private
  public :: read_from_hdf5


  interface read_from_hdf5
     module procedure &
          read_integer, &
!          read_integer_longlong, &
          read_real8, &
          read_integer_array_1d
  end interface read_from_hdf5


contains

  integer(hid_t) function get_native_dtype(ds_id, dname) result(native_dtype)

    integer(hid_t), intent(in) :: ds_id
    character(*),   intent(in) :: dname

    integer(hid_t)             :: dtype_id, native_dtype_id
    integer                    :: class
    integer                    :: ierr
    ! integer :: order, machine_order
    integer(size_t)            :: size_bytes

    !> get the dataset variable type
    !! the "type" and "native_type" are just IDs, the final native type is composed from:
    !! * enddianness
    !! * generic type
    call h5dget_type_f(ds_id, dtype_id, ierr)
    if(ierr/=0) error stop 'h5fortran:reader: get internal dtype ' // dname

    call h5tget_native_type_f(dtype_id, H5T_DIR_ASCEND_F, native_dtype_id, ierr)
    if(ierr/=0) error stop 'h5fortran:reader: get native dtype id ' // dname

    !> we think endianness is handled by HDF5 ... ?
    ! call h5tget_order_f(native_dtype_id, order, ierr)
    ! if(ierr/=0) error stop 'h5fortran:reader: get endianness ' // dname // ' from ' // filename
    ! !> check dataset endianness matches machine (in future, could swap endianness if needed)
    ! call h5tget_order_f(H5T_NATIVE_INTEGER, machine_order, ierr)
    ! if(order /= machine_order) error stop 'h5fortran:reader: endianness does not match machine native ' // dname // ' from ' // filename

    !> compose datatype inferred
    call h5tget_class_f(native_dtype_id, class, ierr)
    if(ierr/=0) error stop 'h5fortran:reader: get class ' // dname

    call h5tget_size_f(native_dtype_id, size_bytes, ierr)
    if(ierr/=0) error stop 'h5fortran:reader: get byte size ' // dname

    call h5tclose_f(dtype_id, ierr)
    if(ierr/=0) error stop 'h5fortran:reader: closing dtype ' // dname

    ! write(*,*) "Class: ", class
    ! write(*,*) "H5T_INTEGER_F: ", H5T_INTEGER_F
    ! write(*,*) "H5T_FLOAT_F:   ", H5T_FLOAT_F
    ! write(*,*) "H5T_STRING_F: ", H5T_STRING_F

    if(class == H5T_INTEGER_F) then
       if(size_bytes == 4) then
          native_dtype = H5T_NATIVE_INTEGER
       elseif(size_bytes == 8) then
          native_dtype = H5T_STD_I64LE
       else
          error stop "h5fortran:reader: expected 32-bit or 64-bit integer:" // dname
       endif
    elseif(class == H5T_FLOAT_F) then
       if(size_bytes == 4) then
          native_dtype = H5T_NATIVE_REAL
       elseif(size_bytes == 8) then
          native_dtype = H5T_NATIVE_DOUBLE
       else
          error stop "h5fortran:reader: expected 32-bit or 64-bit real:" // dname
       endif
    elseif(class == H5T_STRING_F) then
       native_dtype = H5T_NATIVE_CHARACTER
    elseif(class == H5T_ENUM_F) then
       native_dtype = H5T_NATIVE_INTEGER
    else
       error stop "h5fortran:reader: non-handled datatype: " // dname
    endif

  end function get_native_dtype

  subroutine read_integer(x, name, id, got, error)

    integer, intent(out) :: x
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer(hsize_t), parameter  :: xshape(0) = 0
    integer(hid_t) :: dset_id
    integer(hid_t) :: dtype_id

    ! check if dataset exists
    call h5lexists_f(id, name, got, error)
    if (.not. got) then
       error = -1
       write(*,*) 'Dataset ', name, ' does not exist. (?)'
       return
    endif

    ! open dataset
    call h5dopen_f(id, name, dset_id, error)
    if (error /= 0) then
      write(*,'("cannot open hdf5 dataset",/)')
      return
    endif

    ! Get dtype
    dtype_id = get_native_dtype(dset_id, name)



    ! read dataset
    call h5dread_f(dset_id, dtype_id, x, xshape, error)
    if (error /= 0) then
      write(*,'("cannot read hdf5 dataset",/)')
      return
    endif

    ! close dataset
    call h5dclose_f(dset_id, error)
    if (error /= 0) then
      write(*,'("cannot close hdf5 dataset",/)')
      return
    endif

    if (error /= 0) got = .false.
  end subroutine read_integer

  ! subroutine read_integer_longlong(x, name, id, got, error)

  !   integer(kind=8), intent(out) :: x
  !   character(*),    intent(in)  :: name
  !   integer(hid_t),  intent(in)  :: id
  !   logical,         intent(out) :: got
  !   integer,         intent(out) :: error

  !   integer(hsize_t), parameter  :: xshape(0) = 0
  !   integer(hid_t) :: dset_id
  !   integer(hid_t) :: dtype_id

  !   dtype_id = H5T_STD_I64LE

  !   ! check if dataset exists
  !   call h5lexists_f(id, name, got, error)
  !   if (.not.got) return

  !   ! open dataset
  !   call h5dopen_f(id, name, dset_id, error)
  !   if (error /= 0) then
  !     write(*,'("cannot open hdf5 dataset",/)')
  !     return
  !   endif

  !   ! read dataset
  !   call h5dread_f(dset_id, dtype_id, x, xshape, error)
  !   if (error /= 0) then
  !     write(*,'("cannot read hdf5 dataset",/)')
  !     return
  !   endif

  !   ! close dataset
  !   call h5dclose_f(dset_id, error)
  !   if (error /= 0) then
  !     write(*,'("cannot close hdf5 dataset",/)')
  !     return
  !   endif

  !   if (error /= 0) got = .false.
  ! end subroutine read_integer_longlong


  subroutine read_real8(x, name, id, got, error)

    real(kind=8), intent(out) :: x
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer(hsize_t), parameter  :: xshape(0) = 0
    integer(hid_t) :: dset_id
    integer(hid_t) :: dtype_id

    dtype_id = h5t_native_integer

    ! check if dataset exists
    call h5lexists_f(id, name, got, error)
    if (.not.got) return

    ! open dataset
    call h5dopen_f(id, name, dset_id, error)
    if (error /= 0) then
      write(*,'("cannot open hdf5 dataset",/)')
      return
    endif

    ! read dataset
    call h5dread_f(dset_id, dtype_id, x, xshape, error)
    if (error /= 0) then
      write(*,*) '("cannot read hdf5 dataset",/)'
      return
    endif

    ! close dataset
    call h5dclose_f(dset_id, error)
    if (error /= 0) then
      write(*,'("cannot close hdf5 dataset",/)')
      return
    endif

    if (error /= 0) got = .false.
  end subroutine read_real8

  subroutine read_integer_array_1d(x, name, id, got, error)

    integer,         intent(out) :: x(:)
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer, parameter :: ndims = 1
    integer(hsize_t)   :: xshape(ndims)
    integer(hid_t)     :: dset_id
    integer(hid_t)     :: dtype_id

    xshape = shape(x)

    ! check if dataset exists
    call h5lexists_f(id, name, got, error)
    if (.not.got) return

    ! open dataset
    call h5dopen_f(id, name, dset_id, error)
    if (error /= 0) then
      write(*,'("cannot open hdf5 dataset",/)')
      return
    endif

    ! Get dtype
    dtype_id = get_native_dtype(dset_id, name)

    ! read dataset
    call h5dread_f(dset_id, dtype_id, x, xshape, error)
    if (error /= 0) then
      write(*,'("cannot read hdf5 dataset",/)')
      return
    endif

    ! close dataset
    call h5dclose_f(dset_id, error)
    if (error /= 0) then
      write(*,'("cannot close hdf5 dataset",/)')
      return
    endif

    if (error /= 0) got = .false.
  end subroutine read_integer_array_1d

  subroutine read_integer_array_4d(x, name, id, got, error)

    integer,         intent(out) :: x(:,:,:,:)
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer, parameter :: ndims = 4
    integer(hsize_t)   :: xshape(ndims)
    integer(hid_t)     :: dset_id
    integer(hid_t)     :: dtype_id

    xshape = shape(x)

    ! check if dataset exists
    call h5lexists_f(id, name, got, error)
    if (.not.got) return

    ! open dataset
    call h5dopen_f(id, name, dset_id, error)
    if (error /= 0) then
      write(*,'("cannot open hdf5 dataset",/)')
      return
    endif

    ! Get native dtype
    dtype_id = get_native_dtype(dset_id, name)

    ! read dataset
    call h5dread_f(dset_id, dtype_id, x, xshape, error)
    if (error /= 0) then
      write(*,'("cannot read hdf5 dataset",/)')
      return
    endif

    ! close dataset
    call h5dclose_f(dset_id, error)
    if (error /= 0) then
      write(*,'("cannot close hdf5 dataset",/)')
      return
    endif

    if (error /= 0) got = .false.
  end subroutine read_integer_array_4d

end module hdf5_utils
