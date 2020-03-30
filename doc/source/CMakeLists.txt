add_subdirectory(doxygen)

set(DOXYGEN_XML_PATH ${CMAKE_CURRENT_BINARY_DIR}/../build/doxygen/xml)

file(GLOB files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.py")
foreach(file ${files})
  configure_file(${file} ${file} @ONLY)
endforeach()

file(GLOB files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.rst")
foreach(file ${files})
  configure_file(${file} ${file} @ONLY)
endforeach()

