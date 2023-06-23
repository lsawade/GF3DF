module hdf5_utils

  use hdf5
  use iso_c_binding, only:c_loc, c_ptr
  implicit none

  private
  public :: read_from_hdf5, get_dset_dims, get_dset_rank


  interface read_from_hdf5
     module procedure &
          get_dset_rank, &
          get_dset_dims, &
          read_real4_array_2d, &
          read_real8, &
          read_real8_array_1d, &
          read_integer, &
          read_integer_array_1d, &
          read_integer_array_2d, &
          read_integer_array_4d, &
          read_integer_long_array_4d, &
          read_integer_long, &
          read_real_array_2d, &
          read_real_array_5d, &
          read_char_array_1d
  end interface read_from_hdf5


contains

    subroutine get_dset_rank(file_id, name, rank, got, errorflag)

    ! Inputs
    integer(hid_t), intent(in)           :: file_id
    character(*),   intent(in)           :: name

    ! Outputs
    integer, intent(out)                 :: rank
    logical, intent(out)                 :: got
    integer, intent(out)                 :: errorflag

    ! Local
    integer(hid_t)                       :: dset_id, dspace_id

    ! Open dataset
    ! check if dataset exists
    call h5lexists_f(file_id, name, got, errorflag)
    if (.not. got) then
       errorflag = -1
       write(*,*) 'Dataset ', name, ' does not exist. (?)'
       return
    endif

    ! open dataset
    call h5dopen_f(file_id, name, dset_id, errorflag)
    if (errorflag /= 0) then
      write(*,*) "Cannot open hdf5 dataset in get_rank"
      return
    endif

    ! Get the dataspace ID
    call h5dget_space_f(dset_id, dspace_id, errorflag)
    if (errorflag /= 0) then
      write(*,*) "Cannot get dataspace id in get_rank"
      return
    endif


    ! Getting dims from dataspace
    call h5sget_simple_extent_ndims_f(dspace_id, rank, errorflag)
    if (errorflag /= 0) then
      write(*,*) "Cannot get rank in get_rank"
      return
    endif

    ! close dataset
    call h5dclose_f(dset_id, errorflag)
    if (errorflag /= 0) then
      write(*,*) "cannot close hdf5 dataset in get_rank"
      return
    endif

  end subroutine get_dset_rank


  subroutine get_dset_dims(file_id, name, dims, maxdims, got, errorflag)

    ! Inputs
    integer(hid_t), intent(in)           :: file_id
    character(*),   intent(in)           :: name

    ! Outputs
    integer, dimension(:), intent(out)   :: dims, maxdims
    logical, intent(out)                 :: got
    integer, intent(out)                 :: errorflag

    ! Local
    integer(hid_t)                              :: dset_id, dspace_id
    integer(hsize_t), dimension(:), allocatable :: h5dims, h5maxdims

    allocate(h5dims(size(dims)))
    allocate(h5maxdims(size(dims)))

    ! Open dataset
    ! check if dataset exists
    call h5lexists_f(file_id, name, got, errorflag)
    if (.not. got) then
       errorflag = -1
       write(*,*) 'Dataset ', name, ' does not exist. (?)'
       return
    endif

    ! open dataset
    call h5dopen_f(file_id, name, dset_id, errorflag)
    if (errorflag /= 0) then
      write(*,*) "Cannot open hdf5 dataset in get_rank"
      return
    endif

    ! Get the dataspace ID
    call h5dget_space_f(dset_id, dspace_id, errorflag)
    if (errorflag /= 0) then
      write(*,*) "Cannot get dataspace id in get_dims"
      return
    endif


    ! Getting dims from dataspace
    call h5sget_simple_extent_dims_f(dspace_id, h5dims, h5maxdims, errorflag)
    if (errorflag < 0) then
      write(*,*) "Cannot get dims in get_dims"
      write(*,*) "dims", h5dims, "maxdims", h5maxdims
      return
    else
      errorflag = 0
    endif


    ! Allocate memory for the array.
    dims(:) = h5dims(:)
    maxdims(:) = h5maxdims(:)

    ! Deallocate
    deallocate(h5dims)
    deallocate(h5maxdims)

    ! close dataset
    call h5dclose_f(dset_id, errorflag)
    if (errorflag /= 0) then
      write(*,*) "cannot close hdf5 dataset in get_dims"
      return
    endif

  end subroutine get_dset_dims

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

    ! write(*,*) "Class: ", class, size_bytes
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

  subroutine read_integer_long(x, name, id, got, error)

    integer(kind=8), intent(out) :: x
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
  end subroutine read_integer_long

  subroutine read_real4_array_2d(x, name, id, got, error)

    real(kind=4),    intent(out) :: x(:,:)
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer, parameter :: ndims = 2
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
  end subroutine read_real4_array_2d

  subroutine read_real8(x, name, id, got, error)

    real(kind=8),    intent(out) :: x
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer(hsize_t), parameter  :: xshape(0) = 0
    integer(hid_t) :: dset_id
    integer(hid_t) :: dtype_id


    ! check if dataset exists
    call h5lexists_f(id, name, got, error)
    if (.not.got) return

    ! open dataset
    call h5dopen_f(id, name, dset_id, error)
    if (error /= 0) then
      write(*,'("cannot open hdf5 dataset",/)')
      return
    endif

    ! Get datatype
    dtype_id = get_native_dtype(dset_id, name)

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

  subroutine read_real8_array_1d(x, name, id, got, error)

    real(kind=8),    intent(out) :: x(:)
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
  end subroutine read_real8_array_1d

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

  subroutine read_integer_array_2d(x, name, id, got, error)

    integer,    intent(out) :: x(:,:)
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer, parameter :: ndims = 2
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
  end subroutine read_integer_array_2d

  subroutine read_integer_long_array_4d(x, name, id, got, error)

    integer(kind=8),         intent(out)         :: x(:,:,:,:)
    integer,         allocatable, target :: rdata(:,:,:,:,:)
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer, parameter :: ndims = 4
    integer(hsize_t), dimension(ndims)   :: xshape
    integer(hsize_t), dimension(ndims)   :: dims, maxdims
    integer(hsize_t), dimension(1)       :: sdims, smaxdims
    integer(hid_t)     :: dset_id
    integer(hid_t)     :: filetype, memtype, space
    integer(hid_t)     :: dtype_id
    integer            :: i,j

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

    ! Reading the array
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
  end subroutine read_integer_long_array_4d


  subroutine read_integer_array_4d(x, name, id, got, error)

    integer,         intent(out)         :: x(:,:,:,:)
    integer,         allocatable, target :: rdata(:,:,:,:,:)
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer, parameter :: ndims = 4
    integer(hsize_t), dimension(ndims)   :: xshape
    integer(hsize_t), dimension(ndims)   :: dims, maxdims
    integer(hsize_t), dimension(1)       :: sdims, smaxdims
    integer(hid_t)     :: dset_id
    integer(hid_t)     :: filetype, memtype, space
    integer(hid_t)     :: dtype_id
    integer            :: i,j

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

    ! Reading the array
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


    write(*,*) 'post close'

    if (error /= 0) got = .false.
  end subroutine read_integer_array_4d

  subroutine read_real_array_2d(x, name, id, got, error)

    real(kind=8),    intent(out) :: x(:,:)
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer, parameter :: ndims = 2
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
  end subroutine read_real_array_2d

  subroutine read_real_array_5d(x, name, id, got, error)

    real(kind=8),            intent(out) :: x(:,:,:,:,:)
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer, parameter :: ndims = 5
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
  end subroutine read_real_array_5d

  subroutine read_char_array_1d(x, name, id, got, error)

    character(*),intent(out) :: x(:)
    character(*),    intent(in)  :: name
    integer(hid_t),  intent(in)  :: id
    logical,         intent(out) :: got
    integer,         intent(out) :: error

    integer, parameter :: ndims = 1
    integer(hsize_t)   :: xshape(ndims)
    integer(hid_t)     :: dset_id, filetype, memtype
    integer(hid_t)     :: dtype_id, space
    integer(size_t) :: size
    integer(size_t), parameter :: sdim = 32
    integer(hsize_t) :: dims(1)
    integer(hsize_t) :: maxdims(1)
    character(len=sdim), dimension(:), allocatable, target :: rdata(:)
    type(C_PTR) :: f_ptr

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

    ! ! Get dtype
    dtype_id = get_native_dtype(dset_id, name)

    ! Get the datatype and its size.
    call h5dget_type_f(dset_id, filetype, error)

    write(*,*) "filetype", filetype
    call h5tget_size_f(filetype, size, error)
    if (error /= 0) then
      write(*,'("cannot get hdf5 datatype or size",/)')
      return
    endif

    if (size > sdim+1) then
      write(*,*) 'error: character len is too small'
      stop
    endif

    write(*,*) "size", size

    ! Get dataspace.
    call h5dget_space_f(dset_id, space, error)
    if (error /= 0) then
      write(*,'("cannot get hdf5 dataspace",/)')
      return
    endif

    write(*,*) "space", space
    call h5sget_simple_extent_dims_f(space, dims, maxdims, error)
    if ((error < 0)) then
      write(*,'("cannot get hdf5 extent",/)')
      return
    endif

    write (*,*) "dims", dims, "maxdims", maxdims

    allocate(rdata(1:dims(1)))

    ! Create the memory datatype.
    call h5tcopy_f(h5t_fortran_s1,memtype, error)
    call h5tset_size_f(memtype, sdim, error)
    if (error /= 0) then
      write(*,'("cannot get HDF5 memory datatype",/)')
      return
    endif

    ! Read the data.
    f_ptr = C_LOC(rdata(1:dims(1)))

    ! read dataset
    call h5dread_f(dset_id, memtype, f_ptr, error, space)

    ! call h5dread_f(dset_id, memtype, x, xshape, error)
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

    x(:) = rdata(:)

    if (error /= 0) got = .false.
  end subroutine read_char_array_1d

end module hdf5_utils
