add_library(hex_nodes hex_nodes.F90)

target_link_libraries(
    hex_nodes
    constants
)

add_library(reduce reduce.F90)

target_link_libraries(
    reduce
    constants
)


add_library(topo topo.F90)

target_link_libraries(
    topo
    constants
)

add_library(spline_routines spline_routines.F90)

target_link_libraries(
    spline_routines
    constants
)

add_library(rthetaphi_xyz rthetaphi_xyz.F90)

target_link_libraries(
    rthetaphi_xyz
    constants
)


add_library(jacobian jacobian.F90)

target_link_libraries(
    jacobian
    constants
)

add_library(point_location point_location.F90)

target_link_libraries(
    point_location
    ctypes
    constants
    jacobian
)

add_library(gll_library gll_library.F90)

add_library(lagrange_poly lagrange_poly.F90)

target_link_libraries(
    lagrange_poly
    gll_library
)
