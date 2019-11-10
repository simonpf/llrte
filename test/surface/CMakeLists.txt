cmake_minimum_required (VERSION 3.5)

include_directories("${PROJECT_SOURCE_DIR}/src/")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/src/")

add_executable(test_black_surface test_black_surface.cxx)
