function(CUDA_EXECUTABLE)
  set(ONEVALUEARGS NAME FILE)
  set(MULTIVALUEARGS INCLUDES LIBRARIES)
  cmake_parse_arguments(
    CUDA
    ""
    "${ONEVALUEARGS}"
    "${MULTIVALUEARGS}"
    ${ARGN})

  # Include dirs
  list(APPEND CUDA_INCLUDES "-I${CMAKE_SOURCE_DIR}/src")
  list(APPEND CUDA_INCLUDES "-I${CMAKE_SOURCE_DIR}/ext/Catch2/single_include/catch2")

  list(APPEND CUDA_LIBRARIES "-L/usr/local/cuda-10.1/lib64/ -lcudart")

  # Architectures
  set(ARCHS "--cuda-gpu-arch=sm_60")

  set(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${CUDA_NAME}")

  set(FILENAME "${CMAKE_CURRENT_SOURCE_DIR}/${CUDA_FILE}")
  add_custom_target(${CUDA_NAME}
    COMMAND
    clang++ -std=c++2a -DCUDA ${FILENAME} --cuda-path=/usr/local/cuda-10.1 ${CUDA_INCLUDES} ${ARCHS} ${CUDA_LIBRARIES} ${CUDA_LIBRARIES}  -o ${OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Building CUDA executable.")
endfunction(CUDA_EXECUTABLE)
