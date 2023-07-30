# # Library to read/write binary. Main use is to write SAC files.
# add_library(binary_c_io binary_c_io.c ${HEADER_FILES})
# target_include_directories(binary_c_io PUBLIC ${PROJECT_SOURCE_DIR}/include)

# # C library to check file/dir status
# add_library(file_info file_info.c)

# Error function library. Used to make Heaviside functions to convolve with
# the stored traces.

# # Create constants
add_library(constants constants.F90)

# FFT routines from netlib
add_subdirectory(fftpack-sources)

# Utilities of various origin
add_subdirectory(utils)
target_link_libraries(
    utils
    constants
)
target_link_libraries(
    kdtree_search
    heap_sort
)
target_link_libraries(
    stf
    fftpack
    error_functions
)
target_link_libraries(
    hdf5_utils
    hdf5_fortran)


# The code for locating a source in the mesh is here
add_subdirectory(mesh)


add_subdirectory(sources)
target_link_libraries(
    sources
    calendar
    utils
)


add_subdirectory(gf)
target_link_libraries(
    gf
    kdtree_search
    hex_nodes
    hdf5_utils
)

# Files for setup and perform source location in the mesh
add_subdirectory(location)

target_link_libraries(
    source_location
    point_location
    kdtree_search
    rthetaphi_xyz
    reduce
    topo
    spline_routines
    sources
)
target_link_libraries(
    point_location
    jacobian
)

# Interpolation routines
add_subdirectory(interpolation)
target_link_libraries(
    interpolation
    lagrange_poly
)

# SAC writing
add_subdirectory(sac)

# Main function library to loaded by external people
add_subdirectory(gf3d)

target_link_libraries(
    gf3d
    utils
    interpolation
    gf
    sources
    sac
    source_location
    stf
)

