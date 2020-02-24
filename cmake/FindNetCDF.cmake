#[==[
This module uses the nc-config binary to determine location
of include files, libraries and linker flags for 

#]==]

find_program (NCCONFIG nc-config)
if (NCCONFIG)

  execute_process(COMMAND ${NCCONFIG} --libdir OUTPUT_VARIABLE _netcdf_lib)
  execute_process(COMMAND ${NCCONFIG} --includedir OUTPUT_VARIABLE _netcdf_include)
  execute_process(COMMAND ${NCCONFIG} --libs OUTPUT_VARIABLE _netcdf_link)
  execute_process(COMMAND ${NCCONFIG} --version OUTPUT_VARIABLE _netcdf_version)

  # Determine version
  string(REGEX REPLACE "netCDF \([0-9]*.[0-9]*.[0-9]*[\-]*[a-zA-Z0-9]*\)" "\\1"
    _netcdf_version "${_netcdf_version}")
  string(STRIP ${_netcdf_version} _netcdf_version)

  set(NetCDF_FOUND "1")
  set(NetCDF_VERSION "${_netcdf_version}")
  set(NetCDF_VERSION "${_netcdf_version}")
  set(NetCDF_INCLUDE_DIR ${_netcdf_include})
  set(NetCDF_LINKER_FLAGS ${_netcdf_link})

  find_library(NetCDF_LIBRARY netcdf PATHS ${_netcdf_lib})
  message(WARNING "Netcdf: ${NetCDF_LIBRARY}")

  unset(_netcdf_lib)
  unset(_netcdf_include)
  unset(_netcdf_link)
  unset(_netcdf_version)
  unset(_netcdf_version_major)
  unset(_netcdf_version_minor)
  unset(_netcdf_version_patch)
  unset(_netcdf_version_note)

  message(WARNING ${NetCDF_VERSION})
endif (NCCONFIG)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NetCDF
  REQUIRED_VARS NetCDF_LIBRARY NetCDF_INCLUDE_DIR
  VERSION_VAR NetCDF_VERSION)
