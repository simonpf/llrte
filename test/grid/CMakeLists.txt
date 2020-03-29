cmake_minimum_required (VERSION 3.5)

enable_testing()
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/ext/Catch2/single_include/catch2")
add_executable(regular regular.cxx)
