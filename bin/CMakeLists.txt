include_directories("${PROJECT_SOURCE_DIR}/src/")
link_libraries(stdc++fs)
#add_executable(test_module compile_module.cxx)

find_package(NetCDF)
if (NetCDF_FOUND)
  add_executable(exercise_1_a exercise_1_a.cxx)
  target_link_libraries(exercise_1_a "netcdf" "curl")
  target_include_directories(exercise_1_a PUBLIC ${NetCDF_INCLUDE_DIRS})
endif(NetCDF_FOUND)

add_executable(exercise_1_b exercise_1_b.cxx)
add_executable(exercise_2_a exercise_2_a.cxx)
add_executable(exercise_2_b exercise_2_b.cxx)
add_executable(exercise_2_c exercise_2_c.cxx)
add_executable(exercise_3_a exercise_3_a.cxx)
add_executable(exercise_3_b_1 exercise_3_b_1.cxx)
add_executable(exercise_3_b_2 exercise_3_b_2.cxx)
add_executable(exercise_3_b_3 exercise_3_b_3.cxx)
add_executable(exercise_4_a exercise_4_a.cxx)
add_executable(exercise_4_ref exercise_4_ref.cxx)
add_executable(exercise_4_b exercise_4_b.cxx)
add_executable(exercise_5_a exercise_5_a.cxx)
add_executable(exercise_5_b exercise_5_b.cxx)

configure_file(plot_results_1_a.py plot_results_1_a.py COPYONLY)
configure_file(plot_results_1_b.py plot_results_1_b.py COPYONLY)
configure_file(plot_results_2_a.py plot_results_2_a.py COPYONLY)
configure_file(plot_results_2_b.py plot_results_2_b.py COPYONLY)
configure_file(plot_results_2_c.py plot_results_2_c.py COPYONLY)
configure_file(plot_results_3.py plot_results_3.py COPYONLY)
configure_file(plot_results_4.py plot_results_4.py COPYONLY)
configure_file(plot_results_5.py plot_results_5.py COPYONLY)
configure_file(plot_results_5_b_backwards.py plot_results_5_b_backwards.py COPYONLY)
configure_file(plot_results_5_3d.py plot_results_5_3d.py COPYONLY)

find_package(PythonInterp)
if (PYTHONINTERP_FOUND)
  add_custom_target(
    exercise_1
    ${PYTHON_EXECUTABLE} plot_results_1_a.py
    COMMAND ${PYTHON_EXECUTABLE} plot_results_1_b.py
    DEPEND exercise_1_a exercise_1_b
    )
  add_custom_target(
    exercise_2
    ${PYTHON_EXECUTABLE} plot_results_2_a.py
    COMMAND ${PYTHON_EXECUTABLE} plot_results_2_b.py
    COMMAND ${PYTHON_EXECUTABLE} plot_results_2_c.py
    DEPEND exercise_2_a
    )
  add_custom_target(
    exercise_3
    ${PYTHON_EXECUTABLE} plot_results_3.py
    DEPEND exercise_3_b
    )
  add_custom_target(
    exercise_4
    ${PYTHON_EXECUTABLE} plot_results_4.py
    DEPEND exercise_4_b
    )
  add_custom_target(
    exercise_5
    ${PYTHON_EXECUTABLE} plot_results_5.py
    DEPEND exercise_5_a
    )
elseif(PYTHONINTERP_FOUND)
  message(STATUS "No python interpreter found.")
endif(PYTHONINTERP_FOUND)

find_package(NetCDF)
if (NetCDF_FOUND)
  add_executable(exercise_5_b_backwards exercise_5_b_backwards.cxx)
  target_link_libraries(exercise_5_b_backwards ${NetCDF_}"netcdf" "curl")
  target_include_directories(exercise_5_b_backwards PUBLIC ${NetCDF_INCLUDE_DIRS})

  add_executable(exercise_5_3d exercise_5_3d.cxx)
  target_link_libraries(exercise_5_3d ${NetCDF_}"netcdf" "curl")
  target_include_directories(exercise_5_3d PUBLIC ${NetCDF_INCLUDE_DIRS})
endif(NetCDF_FOUND)
