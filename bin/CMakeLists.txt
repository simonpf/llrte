cmake_minimum_required (VERSION 3.5)

include_directories("${PROJECT_SOURCE_DIR}/src/")
link_libraries(stdc++fs)
#add_executable(test_module compile_module.cxx)

add_executable(monte_carlo monte_carlo.cxx)

