# add_library(kdtree_sesarch ${CMAKE_CURRENT_LIST_DIR}/kdtree_search.F90)
# add_library(setup_source_location ${CMAKE_CURRENT_LIST_DIR}/setup_source_location.F90)
add_library(source_location ${CMAKE_CURRENT_LIST_DIR}/source_location.F90)
add_library(point_location ${CMAKE_CURRENT_LIST_DIR}/point_location.F90)

# target_link_libraries(
#     sources
#     point_location
#     source_location
#     kdtree_search
#     # setup_source_location
#     # kdtree_search
# )