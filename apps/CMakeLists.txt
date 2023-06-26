
add_executable(gf3d-print-cmt print_cmt.F90)
target_link_libraries(
    gf3d-print-cmt
    sources
    utils
    calendar
)

add_executable(gf3d-write-sac write_sac_seismograms.F90)
target_link_libraries(
    gf3d-write-sac
    gf3d
)

add_executable(gf3d-write-stf write_source_time_functions.F90)
target_link_libraries(
    gf3d-write-stf
    gf3d
)

# add_executable(readgf hdf5_read.F90)

# #the name of the external library located in /home/path/lib/libdir/lib is    libexternlib.so
# target_link_libraries(
#     readgf
#     gf3d
#     utils
#     ctypes
#     hdf5_utils
#     hdf5_fortran-shared
# )

# add_executable(gf3d-locate-source locate_source.F90)
# target_link_libraries(
#     gf3d-locate-source
#     gf3d
#     utils
#     ctypes
#     setup_source_location
#     hdf5_utils
#     hdf5_fortran-shared
# )

# add_executable(gf3d-interpolate-source interpolate_source.F90)
# target_link_libraries(
#     gf3d-interpolate-source
#     gf3d
#     utils
#     ctypes
#     setup_source_location
#     hdf5_utils
#     hdf5_fortran-shared
# )


# add_executable(test-utils  test_utils.F90)
# target_link_libraries(
#     test-utils
#     utils
# )


# add_executable(test-rw-c test_rw.c)
# target_link_libraries(
#     test-rw-c
#     ${MPI_C_LIBRARIES}
#     hdf5-shared
# )

# add_executable(test-rw-f test_rw.F90)

# target_link_libraries(
#     test-rw-f
#     ${MPI_C_LIBRARIES}
#     hdf5_fortran-shared
# )

