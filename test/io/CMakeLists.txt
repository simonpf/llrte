enable_testing()
include_directories(
  "${PROJECT_SOURCE_DIR}/src/"
  "${PROJECT_SOURCE_DIR}/ext/Catch2/include"
  )

find_package(NetCDF)
if (NetCDF_FOUND)
  add_executable(test_netcdf test_netcdf.cxx)
  target_link_libraries(test_netcdf "netcdf" "curl")
  target_include_directories(test_netcdf PUBLIC ${NetCDF_INCLUDE_DIRS})
endif(NetCDF_FOUND)
