macro (CUDA_EXECUTABLE NAME FILE)

  # Include dirs
  set(INCLUDES "-I${CMAKE_SOURCE_DIR}/src")
  list(APPEND INCLUDES "-I${CMAKE_SOURCE_DIR}/ext/Catch2/single_include/catch2")

  # Architectures
  set(ARCHS "--cuda-gpu-arch=sm_60")

  set(FILENAME "${CMAKE_CURRENT_SOURCE_DIR}/${FILE}")
  add_custom_target(${NAME}
    COMMAND
    clang++ -std=c++2a -DCUDA ${FILENAME} ${INCLUDES} ${ARCHS} -lcudart -o ${NAME}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Building CUDA executable.")
endmacro (CUDA_EXECUTABLE)
