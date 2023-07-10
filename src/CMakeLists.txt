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

# The code for locating a source in the mesh is here
add_subdirectory(mesh)


add_subdirectory(sources)
target_link_libraries(
    sources
    calendar
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
    # sources
    # gf
    # fftpack
    # interpolation
    # sac
    # # location directory
    # kdtree_search
    # point_location
    # source_location
    # # mesh directory
    # gll_library
    # hex_nodes
    # jacobian
    # lagrange_poly
    # reduce
    # rthetaphi_xyz
    # spline_routines
    # topo
    # # Utils directory
    # calendar
    # hdf5_utils
    # heap_sort
    # stf
    # error_functions


# # Create HDF5 utilities
# add_library(hdf5_utils hdf5_utils.F90)
# target_link_libraries(hdf5_utils hdf5_fortran-shared)

# # Create types
# add_library(ctypes ctypes.F90)

# target_link_libraries(
#     ctypes
#     calendar
# )
# # Create utilities
# add_library(utils utils.F90)

# # Link source location
# target_link_libraries(
#     utils
#     ctypes
#     constants
#     hdf5_utils
# )

# # Create sorting algorithm
# add_library(heap_sort heap_sort.F90)

# # Create kdtree search
# add_library(kdtree_search kdtree_search.F90)

# target_link_libraries(
#     kdtree_search
#     heap_sort
# )

# add_subdirectory(io)


# # Create constants
# add_library(setup_source_location setup_source_location.F90)

# # Add CMT library
# add_library(get_cmt get_cmt.F90)

# add_library(source_location source_location.F90)

# # target_link_libraries(
# #     source_location
# #     constants
# #     ctypes
# #     get_cmt
# #     kdtree_search
# #     point_location
# #     rthetaphi_xyz
# #     reduce
# #     topo
# #     utils
# #     spline_routines
# # )

# # Interpolation library
# add_library(interpolation interpolation.F90)

# # target_link_libraries(
# #     interpolation
# #     ctypes
# #     constants
# #     lagrange_poly
# #     utils
# # )

# add_subdirectory(io)

# # Link the libraries.
# # target_link_libraries(
# #     io
# #     utils
# #     ctypes
# #     constants
# #     gll_library
# #     lagrange_poly
# #     setup_source_location
# #     hdf5_utils
# #     hdf5_fortran-shared
# # )

# # Source time function library
# add_library(stf stf.F90)
# target_link_libraries(
#     stf
#     error_functions
# )

# # GF3D is the main module but it will contain some source files
# # To larger functions that are outside the main module
# # otherwise the filee become very bloated.
# set(GF3D_SOURCES ${PROJECT_SOURCE_DIR}/src/gf3d_sub)
# add_library(gf3d gf3d.F90 ${GF3D_SOURCES})

# # Link the libraries.
# target_link_libraries(
#     gf3d
#     io
#     utils
#     ctypes
#     constants
#     setup_source_location
#     source_location
#     interpolation
#     hdf5_utils
#     hdf5_fortran-shared
# )

