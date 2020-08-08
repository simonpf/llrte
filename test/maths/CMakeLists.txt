enable_testing()
include_directories("${PROJECT_SOURCE_DIR}/src/")

get_target_property(Eigen3_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
set(Eigen3_INCLUDE_DIRS ${Eigen3_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS}/unsupported)
include_directories("${PROJECT_SOURCE_DIR}/src/" ${Eigen3_INCLUDE_DIRS})

add_executable(rotations rotations.cxx)
add_executable(geometry geometry.cxx)
