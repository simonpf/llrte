cmake_minimum_required (VERSION 3.5)

include_directories("${PROJECT_SOURCE_DIR}/src/")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/src/")

enable_testing()
add_executable(black_surface black_surface.cxx)
add_test(black_surface black_surface)
add_executable(periodic_boundary periodic_boundary.cxx)
add_test(periodic_boundary periodic_boundary)