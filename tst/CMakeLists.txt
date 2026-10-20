### CMakeLists.txt ---

enable_testing()

find_package(Catch2 3 REQUIRED)

# Shared test helpers (entropy_test_helpers.h, LayerClaim.h) live in tst/igor.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/igor)

add_subdirectory(igor/Core)
add_subdirectory(igor/Math)
add_subdirectory(igor/Alignment)
add_subdirectory(igor/Model)
add_subdirectory(igor/Generation)
add_subdirectory(igor/Inference)
add_subdirectory(igor/Streaming)
add_subdirectory(utils)

######################################################################
### CMakeLists.txt ends here
