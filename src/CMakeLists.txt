# Library to read/write binary. Main use is to write SAC files.
add_library(binary_c_io binary_c_io.c ${HEADER_FILES})
target_include_directories(binary_c_io PUBLIC ${PROJECT_SOURCE_DIR}/include)

# C library to check file/dir status
add_library(file_info file_info.c)

# Error function library. Used to make Heaviside functions to convolve with
# the stored traces.
add_library(error_functions error_functions.F90)

# Create constants
add_library(constants constants.F90)

# Create constants
add_library(calendar calendar.F90)

# The code for locating a source in the mesh is here
add_subdirectory(meshlocation)

# FFT routines from netlib
add_subdirectory(fftpack-sources)

# add_library(fftpack ${FFTPACK_SOURCES})

# Create HDF5 utilities
add_library(hdf5_utils hdf5_utils.F90)
target_link_libraries(hdf5_utils hdf5_fortran-shared)

# Create types
add_library(ctypes ctypes.F90)

target_link_libraries(
    ctypes
    calendar
)
# Create utilities
add_library(utils utils.F90)

# Link source location
target_link_libraries(
    utils
    ctypes
    constants
    hdf5_utils
)

# Create sorting algorithm
add_library(heap_sort heap_sort.F90)

# Create kdtree search
add_library(kdtree_search kdtree_search.F90)

target_link_libraries(
    kdtree_search
    heap_sort
)

# Create constants
add_library(setup_source_location setup_source_location.F90)

# Link source location
target_link_libraries(
    setup_source_location
    ctypes
    kdtree_search
    constants
    hex_nodes
)

# Add CMT library
add_library(get_cmt get_cmt.F90)

add_library(source_location source_location.F90)

target_link_libraries(
    source_location
    constants
    ctypes
    get_cmt
    kdtree_search
    point_location
    rthetaphi_xyz
    reduce
    topo
    utils
    spline_routines
)

# Interpolation library
add_library(interpolation interpolation.F90)

target_link_libraries(
    interpolation
    ctypes
    constants
    lagrange_poly
    utils
)


# Create utilities
add_library(io io.F90 ${HEADER_FILES})

# Link the libraries.
target_link_libraries(
    io
    binary_c_io
    utils
    ctypes
    constants
    gll_library
    lagrange_poly
    setup_source_location
    hdf5_utils
    hdf5_fortran-shared
)

# Source time function library
add_library(stf stf.F90)
target_link_libraries(
    stf
    error_functions
)

add_library(gf3d gf3d.F90)

# Link the libraries.
target_link_libraries(
    gf3d
    io
    utils
    ctypes
    constants
    setup_source_location
    source_location
    interpolation
    hdf5_utils
    hdf5_fortran-shared
)

