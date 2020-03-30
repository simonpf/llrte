add_subdirectory(source)

find_package(Sphinx)

set(SPHINX_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/source)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/build)

add_custom_target(Sphinx
  COMMAND
  ${SPHINX_EXECUTABLE} -b html
  ${SPHINX_SOURCE} ${SPHINX_BUILD}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating documentation with Sphinx")
