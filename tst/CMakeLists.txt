### CMakeLists.txt ---

enable_testing()

find_package(Catch2 3 REQUIRED)

if(WIN32)
    set(PIXI_BIN_DIR "$ENV{CONDA_PREFIX}/Library/bin")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/igor)

add_subdirectory(igor/Core)
add_subdirectory(igor/Streaming)
add_subdirectory(igor/Math)
add_subdirectory(igor/Model)

######################################################################
### CMakeLists.txt ends here
