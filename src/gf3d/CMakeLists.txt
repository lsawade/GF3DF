add_library(gf3d SHARED
  ${CMAKE_CURRENT_LIST_DIR}/gf3d.F90
  ${CMAKE_CURRENT_LIST_DIR}/seismograms.F90
  ${CMAKE_CURRENT_LIST_DIR}/seismograms_windows.F90
  ${CMAKE_CURRENT_LIST_DIR}/interpolate_source.F90
)
