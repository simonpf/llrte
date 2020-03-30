macro (CUDA_EXECUTABLE NAME FILE NETCDF)

  # Include dirs
  set(INCLUDES "-I${CMAKE_SOURCE_DIR}/src")
  list(APPEND INCLUDES "-I${CMAKE_SOURCE_DIR}/ext/Catch2/single_include/catch2")

  set(LIBRARIES "-lcudart")

  if (${NETCDF})
    find_package(NetCDF)
    list(APPEND INCLUDES "-I${NetCDF_INCLUDE_DIRS}")
    list(APPEND LIBRARIES "-lnetcdf -lcurl")
  endif()

  # Architectures
  set(ARCHS "--cuda-gpu-arch=sm_60")

  set(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${NAME}")

  set(FILENAME "${CMAKE_CURRENT_SOURCE_DIR}/${FILE}")
  add_custom_target(${NAME}
    COMMAND
    clang++ -std=c++2a -DCUDA ${FILENAME} ${INCLUDES} ${ARCHS} ${LIBRARIES} -o ${OUTPUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Building CUDA executable.")
endmacro (CUDA_EXECUTABLE)
