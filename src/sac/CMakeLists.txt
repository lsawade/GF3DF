# Binary IO library for writing SAC
add_library(binary_c_io binary_c_io.c ${HEADER_FILES})
target_include_directories(binary_c_io PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Add sac writing library
add_library(sac ${CMAKE_CURRENT_LIST_DIR}/sac.F90)
target_link_libraries(sac binary_c_io)
