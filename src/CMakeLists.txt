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

# Create utilities
add_library(io io.F90)

# Link the libraries.
target_link_libraries(
    io
    utils
    ctypes
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
    hdf5_utils
    hdf5_fortran-shared
)

