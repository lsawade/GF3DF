add_library(interpolate_seismograms ${CMAKE_CURRENT_LIST_DIR}/interpolate_seismograms.F90)

add_library(seismograms ${CMAKE_CURRENT_LIST_DIR}/seismograms.F90)
target_link_libraries(
  seismograms
  interpolate_seismograms
)

add_library(seismograms_windows ${CMAKE_CURRENT_LIST_DIR}/seismograms_windows.F90)

add_library(gf3d SHARED ${CMAKE_CURRENT_LIST_DIR}/gf3d.F90)

target_link_libraries(
  gf3d
  seismograms
  seismograms_windows
)