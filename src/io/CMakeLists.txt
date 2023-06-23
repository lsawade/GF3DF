# Library to read/write binary. Main use is to write SAC files.
add_library(binary_c_io binary_c_io.c ${HEADER_FILES})
target_include_directories(binary_c_io PUBLIC ${PROJECT_SOURCE_DIR}/include)

# C library to check file/dir status
add_library(file_info file_info.c)

add_library(io "")

target_sources(io
  INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}/GF.f90
    ${CMAKE_CURRENT_LIST_DIR}/cmt.f90
    ${CMAKE_CURRENT_LIST_DIR}/output.f90
    ${CMAKE_CURRENT_LIST_DIR}/sac.f90
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/io.F90
  )

# target_link_libraries(
#     io
#     binary_c_io
# )

# # Create utilities
# add_library(io "")

# # Link the libraries.
# target_link_libraries(
#     io
#     binary_c_io
#     utils
#     ctypes
#     constants
#     gll_library
#     lagrange_poly
#     setup_source_location
#     hdf5_utils
#     hdf5_fortran-shared
# )