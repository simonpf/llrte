cmake_minimum_required (VERSION 3.5)
add_subdirectory(source)

find_package(Sphinx REQUIRED)

set(SPHINX_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/source)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/build)

set(DOXYGEN_XML_PATH ${CMAKE_CURRENT_BINARY_DIR/build/doxygen/xml})

file(GLOB files "*.py" LIST_DIRECTORIES false)
foreach(file ${files})
  configure_file(${i}, ${i} @ONLY)
endforeach()

file(GLOB files ".rst" LIST_DIRECTORIES false)
foreach(file ${files})
  configure_file(${i}, ${i} @ONLY)
endforeach()

add_custom_target(Sphinx ALL
  COMMAND
  ${SPHINX_EXECUTABLE} -b html
  ${SPHINX_SOURCE} ${SPHINX_BUILD}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating documentation with Sphinx")
