add_library(gf
  ${CMAKE_CURRENT_LIST_DIR}/gf.F90
  ${CMAKE_CURRENT_LIST_DIR}/io.F90
  # ${CMAKE_CURRENT_LIST_DIR}/types.F90
)

target_link_libraries(
  gf
  kdtree_search
)
