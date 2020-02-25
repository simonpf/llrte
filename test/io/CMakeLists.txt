cmake_minimum_required (VERSION 3.5)

enable_testing()
include_directories(
  "${PROJECT_SOURCE_DIR}/src/"
  "${PROJECT_SOURCE_DIR}/ext/"
  )

find_package(NetCDF)
if (NetCDF_FOUND)
  #link_directories(${NetCDF_LIBRARY})
  add_executable(bla test_netcdf.cxx)
  #target_link_options(bla PUBLIC -L/usr/local/lib -lnetcdf -lcurl)
  target_link_libraries(bla ${NetCDF_}"netcdf" "curl")
  target_include_directories(bla PUBLIC ${NetCDF_INCLUDE_DIRS})
  message(WARNING ${NetCDF_LIBRARY})
endif(NetCDF_FOUND)
