add_library(utils           ${CMAKE_CURRENT_LIST_DIR}/utils.F90)
add_library(calendar        ${CMAKE_CURRENT_LIST_DIR}/calendar.F90)
add_library(hdf5_utils      ${CMAKE_CURRENT_LIST_DIR}/hdf5_utils.F90)
add_library(heap_sort       ${CMAKE_CURRENT_LIST_DIR}/heap_sort.F90)
add_library(stf             ${CMAKE_CURRENT_LIST_DIR}/stf.F90)
add_library(error_functions ${CMAKE_CURRENT_LIST_DIR}/error_functions.F90)

add_library(kdtree_search
  ${CMAKE_CURRENT_LIST_DIR}/kdtree_search.F90
)
