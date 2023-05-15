cmake_minimum_required (VERSION 3.6)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

project (GF3D
    VERSION 0.1
    DESCRIPTION "Library to extract Green functions")

enable_language (C Fortran)

# set (FFLAGS "-Og -g -Wall -Wextra -pedantic -std=f95 -fimplicit=none -fcheck-all -fbacktrace")

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Fix where you are
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

# Set library output
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/archive)

message(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}")

# Set binary output
set (CMAKE_BINARY_DIR ${PROJECT_BINARY_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# For this make sure that the environment variable HDF5_ROOT is set!
find_package(HDF5 REQUIRED)
find_package(HDF5 COMPONENTS C Fortran HL)


# Set include directories
include_directories(
    ${HDF5_INCLUDE_DIRS}
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

message(HDF5_INCLUDE_DIRS="${HDF5_INCLUDE_DIRS}")
message(HDF5_LIBRARIES="${HDF5_LIBRARIES}")


# The compiled library code is here
add_subdirectory(src)

# The executable code is here
add_subdirectory(apps)





