get_target_property(Eigen3_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
set(Eigen3_INCLUDE_DIRS ${Eigen3_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS}/unsupported)
include_directories("${PROJECT_SOURCE_DIR}/src/" ${Eigen3_INCLUDE_DIRS})

configure_file(data/scattering_data_azimuthally_random.nc
  data/scattering_data_azimuthally_random.nc COPYONLY)
configure_file(test_point_scattering.py test_point_scattering.py COPYONLY)

add_pxx_module(
  SOURCE gridded_optical_properties.h
  MODULE gridded_optical_properties
  INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/src ${Eigen3_INCLUDE_DIRS}
  )
