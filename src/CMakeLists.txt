# Create constants
add_library(constants constants.F90)


# The compiled library code is here
add_subdirectory(meshlocation)

# Create HDF5 utilities
add_library(hdf5_utils hdf5_utils.F90)
target_link_libraries(hdf5_utils hdf5_fortran-shared)

# Create types
add_library(ctypes ctypes.F90)

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
add_library(io io.F90)

# Link the libraries.
target_link_libraries(
    io
    utils
    ctypes
    constants
    gll_library
    lagrange_poly
    hdf5_utils
    hdf5_fortran-shared
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

