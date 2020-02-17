cmake_minimum_required (VERSION 3.5)

enable_testing()
include_directories("${PROJECT_SOURCE_DIR}/src/")

find_package(NetCDF)
if (NetCDF_FOUND)
  add_executable(test_netcdf test_netcdf.cxx)
  target_include_directories(test_netcdf PRIVATE ${NetCDF_INCLUDE_DIRS})
  target_link_libraries(test_netcdf ${NetCDF_LIBRARIES})
endif(NetCDF_FOUND)
