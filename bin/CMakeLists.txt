cmake_minimum_required (VERSION 3.5)

include_directories("${PROJECT_SOURCE_DIR}/src/")
link_libraries(stdc++fs)
#add_executable(test_module compile_module.cxx)

add_executable(exercise_1_a exercise_1_a.cxx)
add_executable(exercise_1_b exercise_1_b.cxx)
add_executable(exercise_2_a exercise_2_a.cxx)
add_executable(exercise_2_b exercise_2_b.cxx)
add_executable(exercise_3_a exercise_3_a.cxx)

configure_file(plot_results_1_a.py plot_results_1_a.py COPYONLY)
configure_file(plot_results_1_b.py plot_results_1_b.py COPYONLY)
configure_file(plot_results_2_a.py plot_results_2_a.py COPYONLY)
configure_file(plot_results_2_b.py plot_results_2_b.py COPYONLY)
configure_file(plot_results_3_a.py plot_results_3_a.py COPYONLY)

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
    DEPEND exercise_2_a
    )
elseif(PYTHONINTERP_FOUND)
  message(STATUS "No python interpreter found.")
endif(PYTHONINTERP_FOUND)
