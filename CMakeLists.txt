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

if(POLICY CMP0079 )
  cmake_policy(SET CMP0079 NEW)
endif()

CMAKE_POLICY(SET CMP0054 NEW) #-- fixes IF() with quoted args
CMAKE_POLICY(SET CMP0057 NEW) #-- enables IF(.. IN_LIST ..)
CMAKE_POLICY(SET CMP0074 NEW) #-- Removes issue about setting HDF5_ROOT

# set (FFLAGS "-fPIC")
# set (FFLAGS "-Og -g -Wall -Wextra -pedantic -std=f95 -fimplicit=none -fcheck-all -fbacktrace")

# Fix where you are
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

# Set library output
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/archive)

message(STATUS "CMAKE_Fortran_MODULE_DIRECTORY: ${CMAKE_Fortran_MODULE_DIRECTORY}")

# Set binary output
set (CMAKE_BINARY_DIR ${PROJECT_BINARY_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

message(STATUS "hello $ENV{HOSTNAME}")

# For this make sure that the environment variable HDF5_ROOT is set!
if(DEFINED ENV{HDF5_ROOT})

    set(HDF5_ROOT "$ENV{HDF5_ROOT}")

    if((DEFINED ENV{HDF5_LIBRARIES}) AND (DEFINED ENV{HDF5_INCLUDE_DIRS}))
        set(HDF5_INCLUDE_DIRS "$ENV{HDF5_INCLUDE_DIRS}")
        set(HDF5_LIBRARIES "$ENV{HDF5_LIBRARIES}")

    # If only one of the HDF5 flags is set, make sure that they are found
    elseif((DEFINED ENV{HDF5_LIBRARIES}) AND (NOT DEFINED ENV{HDF5_INCLUDE_DIRS}))
        message(FATAL_ERROR "Forgot to set HDF5_INCLUDE_DIRS?")
        return()
    elseif((NOT DEFINED ENV{HDF5_LIBRARIES}) AND (DEFINED ENV{HDF5_INCLUDE_DIRS}))
        message(FATAL_ERROR "Forgot to set HDF5_LIBRARIES?")
        return()
    endif()

endif()

find_package(HDF5 REQUIRED COMPONENTS C Fortran HL)


message(STATUS =============== HDF5 setup ===================)
message(STATUS "HDF5_ROOT = ${HDF5_ROOT}")
message(STATUS "HDF5_INCLUDE_DIRS = ${HDF5_INCLUDE_DIRS}")
message(STATUS "HDF5_LIBRARIES = ${HDF5_LIBRARIES}")
message(STATUS ==============================================)

# Set location of C header files
set(HEADER_FILES ${PROJECT_SOURCE_DIR}/include/config.h)

# if (HDF5_FOUND)

# Set include directories
include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Fix for missing library link
list(GET HDF5_INCLUDE_DIRS 0 HDF5_ROOT)
get_filename_component(HDF5_ROOT "${HDF5_ROOT}" DIRECTORY)
set(HDF5_LIBDIR "${HDF5_ROOT}/lib")

message(STATUS "===========================")
message(STATUS "HDF5_LIBDIR is now: ${HDF5_LIBDIR}")
message(STATUS "===========================")

link_directories(${HDF5_LIBDIR})

# include(FetchContent) # once in the project to include the module

# FetchContent_Declare(fftpack
#                      GIT_REPOSITORY https://github.com/fortran-lang/fftpack
#                      GIT_TAG        113982e)
# FetchContent_MakeAvailable(fftpack)


# The compiled library code is here
add_subdirectory(src)

# The executable code is here
add_subdirectory(apps)





