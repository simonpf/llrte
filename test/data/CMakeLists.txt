cmake_minimum_required (VERSION 3.5)

enable_testing()
include_directories("${PROJECT_SOURCE_DIR}/ext/Catch2/single_include/catch2")
include_directories("${PROJECT_SOURCE_DIR}/src/")
add_executable(test_data test_data.cxx)
add_executable(data data.cxx)
