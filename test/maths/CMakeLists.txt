cmake_minimum_required (VERSION 3.5)

enable_testing()
include_directories("${PROJECT_SOURCE_DIR}/src/")
add_executable(rotations rotations.cxx)
add_executable(geometry geometry.cxx)
