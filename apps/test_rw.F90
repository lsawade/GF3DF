! ************************************************************
!
!  this example shows how to read and write data to a
!  dataset.  the program first writes integers to a dataset
!  with dataspace dimensions of dim0xdim1, then closes the
!  file.  next, it reopens the file, reads back the data, and
!  outputs it to the screen.
!
!  this file is intended for use with hdf5 library verion 1.8
!
! ************************************************************

program main

  use hdf5

  implicit none

!   character(len=14), parameter :: filename = "h5ex_d_rdwr.h5"
!   character(len=3) , parameter :: dataset = "DS1"
  character(len=17), parameter :: filename = "single_element.h5"
  character(len=5) , parameter :: dataset = "ibool"
  integer          , parameter :: dim0     = 1
  integer          , parameter :: dim1     = 3
  integer          , parameter :: dim2     = 3
  integer          , parameter :: dim3     = 3

  integer :: hdferr
  integer(hid_t) :: file, space, dset ! handles
  integer(hsize_t), dimension(1:4)           :: dims = (/dim0, dim1, dim2, dim3/) ! size read/write buffer
  integer(kind=8)         , dimension(1:dim0,1:dim1, 1:dim2, 1:dim3) :: wdata, &  ! write buffer
                                                                rdata     ! read buffer
  integer :: i, j, k, l
  !
  ! initialize fortran interface.
  !
  call h5open_f(hdferr)
  !
  ! initialize data.
  !
!   do i = 1, dim0
!      do j = 1, dim1
!         wdata(i,j) = (i-1)*(j-1)-(j-1)
!      enddo
!   enddo
  !
  ! create a new file using the default properties.
  !

!   call h5fcreate_f(filename, h5f_acc_trunc_f, file, hdferr)
!   !
!   ! create dataspace.  setting size to be the current size.
!   !
!   call h5screate_simple_f(2, dims, space, hdferr)
!   !
!   ! create the dataset.  we will use all default properties for this
!   ! example.
!   !
!   call h5dcreate_f(file, dataset, h5t_std_i32le, space, dset, hdferr)
!   !
!   ! write the data to the dataset.
!   !
!   call h5dwrite_f(dset, h5t_native_integer, wdata, dims, hdferr)
!   !
!   ! close and release resources.
!   !
!   call h5dclose_f(dset , hdferr)
!   call h5sclose_f(space, hdferr)
!   call h5fclose_f(file , hdferr)
  !
  ! now we begin the read section of this example.
  !
  !
  ! open file and dataset using the default properties.
  !
  call h5fopen_f(filename, h5f_acc_rdonly_f, file, hdferr)
  call h5dopen_f (file, dataset, dset, hdferr)
  !
  ! read the data using the default properties.
  !
  call h5dread_f(dset, H5T_STD_I64LE, rdata, dims, hdferr)
  !
  ! output the data to the screen.
  !
  write(*, '(/,a,":")') dataset
  do i=1, dim3
    write(*,'(" [")')
     do j=1, dim2
        write(*,'("   [")', advance='no')
        write(*,'(80i3)', advance='no') rdata(1,:,j,i)
        write(*,'("   ]")')
     enddo
    write(*,'(" ]")')
  enddo
  write(*, '(/)')
  !
  ! close and release resources.
  !
  call h5dclose_f(dset , hdferr)
  call h5fclose_f(file , hdferr)

end program main