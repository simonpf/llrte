enable_testing()
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/ext/Catch2/single_include/catch2")
add_executable(regular regular.cxx)

find_package(NetCDF)
if (NetCDF_FOUND)
  add_executable(interpolation interpolation.cxx)
  target_link_libraries(interpolation "netcdf" "curl")
  target_include_directories(interpolation PUBLIC ${NetCDF_INCLUDE_DIRS})
endif(NetCDF_FOUND)

configure_file(plot_results_interpolation.py plot_results_interpolation.py COPYONLY)
