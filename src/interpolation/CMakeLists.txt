add_library(interpolation
  ${CMAKE_CURRENT_LIST_DIR}/interpolation.F90
)
add_library(splines ${CMAKE_CURRENT_LIST_DIR}/splines.F90)
target_link_libraries(
  interpolation
  splines
)