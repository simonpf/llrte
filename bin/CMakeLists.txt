cmake_minimum_required (VERSION 3.5)

include_directories("${PROJECT_SOURCE_DIR}/src/")
link_libraries(stdc++fs)
#add_executable(test_module compile_module.cxx)

add_executable(exercise_1_a exercise_1_a.cxx)
add_executable(exercise_1_b exercise_1_b.cxx)

configure_file(plot_results_1_a.py plot_results_1_a.py COPYONLY)
configure_file(plot_results_1_b.py plot_results_1_b.py COPYONLY)

find_package(PythonInterp)
if (PYTHONINTERP_FOUND)
  add_custom_target(
    run
    ${PYTHON_EXECUTABLE} plot_results_1_a.py
    COMMAND ${PYTHON_EXECUTABLE} plot_results_1_b.py
    DEPEND exercise_1_a exercise_1_b
    )
elseif(PYTHONINTERP_FOUND)
  message(STATUS "No python interpreter found.")
endif(PYTHONINTERP_FOUND)
