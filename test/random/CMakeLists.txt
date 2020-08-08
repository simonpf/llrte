enable_testing()
include_directories(
  "${PROJECT_SOURCE_DIR}/src/"
  "${PROJECT_SOURCE_DIR}/ext/Catch2/single_include/catch2"
  )

get_target_property(Eigen3_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
set(Eigen3_INCLUDE_DIRS ${Eigen3_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS}/unsupported)
include_directories("${PROJECT_SOURCE_DIR}/src/" ${Eigen3_INCLUDE_DIRS})

find_package(OpenMP)
find_package(NetCDF)
if (NetCDF_FOUND)

  if (OpenMP_FOUND)
  add_executable(test_random random.cxx)
  target_link_libraries(test_random "netcdf" "curl")
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif (OpenMP_FOUND)

  cuda_executable(
    NAME test_random_cuda
    FILE random.cu
    INCLUDES ${NetCDF_INCLUDE_DIRS}
    LIBRARIES -lnetcdf -lcurl)

endif(NetCDF_FOUND)
