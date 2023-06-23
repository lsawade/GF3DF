program test_utils

  use utils, only: throwerror, file_exists, dir_exists, dir_writeable, &
                   delete_file, mkdir

  character(len=250) :: dirname = 'src/'
  character(len=250) :: dirname_bad = 'src_bad'
  character(len=250) :: filename = 'CMTSOLUTION'
  character(len=250) :: filename_bad = 'CMTSOLUTION_bad'
  logical :: gooddirres, bad_dirres

  ! --- Check if directory exists ---
  write (*,*) "DIR  Good", dir_exists(dirname), "Bad", dir_exists(dirname_bad)
  write (*,*) "FILE Good", file_exists(filename), "Bad", file_exists(filename_bad)

  ! --- Check if directory is writable ---
  write (*,*) "DIR WRITEABLE: Good", gooddirres,  "Bad", bad_dirres

  ! --- Check whether I can make directories
  write(*,*) 'Making the bad directory ...'
  call mkdir(dirname_bad)
  write (*,*) "FILE Good", file_exists(filename), "Bad", file_exists(filename_bad)
  write (*,*) "DIR WRITEABLE: Good", dir_writeable(dirname),  "Bad", dir_writeable(dirname_bad)

  write(*,*) 'Deleting the bad directory ...'
  call delete_file(dirname_bad)
  write (*,*) "FILE Good", file_exists(filename), "Bad", file_exists(filename_bad)
  write (*,*) "DIR WRITEABLE: Good", dir_writeable(dirname),  "Bad", dir_writeable(dirname_bad)

end program test_utils
