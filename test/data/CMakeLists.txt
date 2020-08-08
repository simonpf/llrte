enable_testing()
include_directories("${PROJECT_SOURCE_DIR}/ext/Catch2/single_include/catch2")
include_directories("${PROJECT_SOURCE_DIR}/src/")

get_target_property(Eigen3_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
set(Eigen3_INCLUDE_DIRS ${Eigen3_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS}/unsupported)
include_directories("${PROJECT_SOURCE_DIR}/src/" ${Eigen3_INCLUDE_DIRS})

add_executable(test_data test_data.cxx)
add_executable(data data.cxx)

# CUDA test
cuda_executable(
  NAME data_cuda
  FILE data.cu)
