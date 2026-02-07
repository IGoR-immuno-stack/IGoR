enable_testing()
# CMake configuration for tests

# Find Catch2
find_package(Catch2 3 REQUIRED)

# Add test sources
file(GLOB TEST_SOURCES "*.cpp")

add_executable(igor_tests ${TEST_SOURCES})

# Link against the core library and Catch2
# (Assuming igor::Core is the main library target)
target_link_libraries(igor_tests PRIVATE igor::Core Catch2::Catch2WithMain)

# Set include directories if needed
# target_include_directories(igor_tests PRIVATE ../src/igor/Core)

# Register the tests with CTest
include(Catch)
catch_discover_tests(igor_tests)

# --- Math Tests ---
file(GLOB MATH_TEST_SOURCES "igor/Math/*.cpp")

if(MATH_TEST_SOURCES)
    add_executable(math_tests ${MATH_TEST_SOURCES})
    set_target_properties(math_tests PROPERTIES CXX_STANDARD 23)
    target_link_libraries(math_tests PRIVATE igor::Math Catch2::Catch2WithMain)
    include(Catch)
    catch_discover_tests(math_tests)
endif()
