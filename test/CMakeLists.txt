add_subdirectory(io)
add_subdirectory(data)
add_subdirectory(surface)
add_subdirectory(maths)
add_subdirectory(random)
add_subdirectory(grid)
add_subdirectory(scattering)

file(GLOB_RECURSE PYTHON_FILES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
foreach(file ${PYTHON_FILES})
  configure_file(${file} ${file} @ONLY)
endforeach()
