enable_testing()
include_directories(
  "${PROJECT_SOURCE_DIR}/src/"
  "${PROJECT_SOURCE_DIR}/ext/Catch2/include"
  )

find_package(NetCDF)
if (NetCDF_FOUND)
  add_executable(test_random test_random.cxx)
  target_link_libraries(test_random ${NetCDF_}"netcdf" "curl")
  target_include_directories(test_random PUBLIC ${NetCDF_INCLUDE_DIRS})
endif(NetCDF_FOUND)
