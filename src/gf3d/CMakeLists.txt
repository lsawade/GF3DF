add_library(interpolate_seismograms ${CMAKE_CURRENT_LIST_DIR}/interpolate_seismograms.F90)
target_link_libraries(
  interpolate_seismograms
  sources
  utils
)


add_library(seismograms ${CMAKE_CURRENT_LIST_DIR}/seismograms.F90)
target_link_libraries(
  seismograms
  interpolate_seismograms
  sources
  utils
)

add_library(seismograms_windows ${CMAKE_CURRENT_LIST_DIR}/seismograms_windows.F90)
target_link_libraries(
  seismograms_windows
  interpolate_seismograms
  sources
  utils
)

add_library(gf3d SHARED ${CMAKE_CURRENT_LIST_DIR}/gf3d.F90)

target_link_libraries(
  gf3d
  seismograms
  seismograms_windows
)