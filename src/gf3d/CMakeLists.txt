add_library(gf3d SHARED
  ${CMAKE_CURRENT_LIST_DIR}/gf3d.F90
  ${CMAKE_CURRENT_LIST_DIR}/gf3d_get_seismograms.F90
  ${CMAKE_CURRENT_LIST_DIR}/gf3d_write_seismograms.F90
)

target_link_libraries(
  gf3d
  point_location
  sources
)