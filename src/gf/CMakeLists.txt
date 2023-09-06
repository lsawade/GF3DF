# add_library(gfio )

add_library(gf ${CMAKE_CURRENT_LIST_DIR}/gf.F90)

target_link_libraries(
  gf
  kdtree_search
)
