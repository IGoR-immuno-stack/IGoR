### CMakeLists.txt ---

enable_testing()

find_package(Catch2 3 REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/igor)

add_subdirectory(igor/Core)
add_subdirectory(igor/Streaming)

######################################################################
### CMakeLists.txt ends here
