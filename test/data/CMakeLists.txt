cmake_minimum_required (VERSION 3.5)

enable_testing()
include_directories("${PROJECT_SOURCE_DIR}/src/")
add_executable(data test_data.cxx)
