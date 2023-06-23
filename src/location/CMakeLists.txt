add_library(source_location ${CMAKE_CURRENT_LIST_DIR}/source_location.F90)
add_library(setup_source_location ${CMAKE_CURRENT_LIST_DIR}/setup_source_location.F90)
add_library(point_location ${CMAKE_CURRENT_LIST_DIR}/point_location.F90)
add_library(kdtree_search ${CMAKE_CURRENT_LIST_DIR}/kdtree_search.F90)

target_link_libraries(
    source_location
    point_location
    kdtree_search
    sources
    setup_source_location
)