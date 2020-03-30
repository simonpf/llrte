enable_testing()
include_directories(
  "${PROJECT_SOURCE_DIR}/src/"
  "${PROJECT_SOURCE_DIR}/ext/Catch2/single_include/catch2"
  )

find_package(OpenMP)
find_package(NetCDF)
if (NetCDF_FOUND)

  add_executable(test_random test_random.cxx)
  target_link_libraries(test_random ${NetCDF_}"netcdf" "curl")
  target_include_directories(test_random PUBLIC ${NetCDF_INCLUDE_DIRS})

  if (OpenMP_FOUND)
  add_executable(random random.cxx)
  target_link_libraries(random ${NetCDF_}"netcdf" "curl")
  SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif (OpenMP_FOUND)

  cuda_executable(random_cuda random.cu TRUE)

endif(NetCDF_FOUND)
