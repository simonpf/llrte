include_directories("${PROJECT_SOURCE_DIR}/src/")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/src/")

get_target_property(Eigen3_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
set(Eigen3_INCLUDE_DIRS ${Eigen3_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS}/unsupported)
include_directories("${PROJECT_SOURCE_DIR}/src/" ${Eigen3_INCLUDE_DIRS})

enable_testing()
configure_file(plot_results_reflecting.py plot_results_reflecting.py COPYONLY)

find_package(NetCDF)
if (NetCDF_FOUND)
  add_executable(black_surface black_surface.cxx)
  add_test(black_surface black_surface)
  target_link_libraries(black_surface "netcdf" "curl")
  target_include_directories(black_surface PUBLIC ${NetCDF_INCLUDE_DIRS})

  add_executable(reflecting_surface reflecting_surface.cxx)
  add_test(reflecting_surface reflecting_surface)
  target_link_libraries(reflecting_surface "netcdf" "curl")
  target_include_directories(reflecting_surface PUBLIC ${NetCDF_INCLUDE_DIRS})

  add_executable(periodic_boundary periodic_boundary.cxx)
  add_test(periodic_boundary periodic_boundary)
  target_link_libraries(periodic_boundary "netcdf" "curl")
  target_include_directories(periodic_boundary PUBLIC ${NetCDF_INCLUDE_DIRS})
endif(NetCDF_FOUND)
