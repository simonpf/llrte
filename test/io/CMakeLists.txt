enable_testing()
include_directories(
  "${PROJECT_SOURCE_DIR}/src/"
  "${PROJECT_SOURCE_DIR}/ext/Catch2/include"
  )

get_target_property(Eigen3_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
set(Eigen3_INCLUDE_DIRS ${Eigen3_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS}/unsupported)
include_directories("${PROJECT_SOURCE_DIR}/src/" ${Eigen3_INCLUDE_DIRS})

find_package(NetCDF)
if (NetCDF_FOUND)
  add_executable(test_netcdf test_netcdf.cxx)
  target_link_libraries(test_netcdf "netcdf" "curl")
  target_include_directories(test_netcdf PUBLIC ${NetCDF_INCLUDE_DIRS})
  set_property(TARGET test_netcdf PROPERTY CXX_STANDARD 20)
endif(NetCDF_FOUND)

