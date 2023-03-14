cmake_minimum_required (VERSION 3.6)

project (GF3D Fortran)

enable_language (C Fortran)

# Fix where you are
set(CMAKE_Fortran_MODULE_DIRECTORY ./include)

# Set library output
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ./lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ./archive)

# Set binary output
set(BINARY_DIR ./bin)

# For this make sure that the environment variable HDF5_ROOT is set!
find_package(HDF5 REQUIRED)
find_package(HDF5 COMPONENTS Fortran HL)

# Set library directories
link_directories (
    ${HDF5_LIBRARIES}
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
)

# Set include directories
include_directories(
    ${HDF5_INCLUDE_DIRS}
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Main source file
set(main-source_file hdf5_read.F90)

add_library(
    hdf5_utils
    src/hdf5_utils.F90
)

target_link_libraries(
    hdf5_utils
    hdf5_fortran
)

add_executable(gf3d ${main-source_file})


#the name of the external library located in /home/path/lib/libdir/lib is    libexternlib.so
target_link_libraries(
    gf3d
    hdf5_utils
    hdf5_fortran
)