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

# Provide the project source directory so tests can locate data files
target_compile_definitions(igor_tests PRIVATE IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# Set include directories if needed
# target_include_directories(igor_tests PRIVATE ../src/igor/Core)

# Register the tests with CTest
include(Catch)
catch_discover_tests(igor_tests
    PROPERTIES TIMEOUT 3600  # generation test with 1M sequences needs headroom
)

# --- Math Tests ---
file(GLOB MATH_TEST_SOURCES "igor/Math/*.cpp")

if(MATH_TEST_SOURCES)
    add_executable(math_tests ${MATH_TEST_SOURCES})
    set_target_properties(math_tests PROPERTIES CXX_STANDARD 23)
    target_link_libraries(math_tests PRIVATE igor::Math Catch2::Catch2WithMain)
    include(Catch)
    catch_discover_tests(math_tests)
endif()

# --- Model Tests ---
file(GLOB MODEL_TEST_SOURCES "igor/Model/*.cpp")
list(REMOVE_ITEM MODEL_TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/igor/Model/test_SamplingEngine.cpp")

if(MODEL_TEST_SOURCES)
    add_executable(model_tests ${MODEL_TEST_SOURCES})
    set_target_properties(model_tests PROPERTIES CXX_STANDARD 23)
    target_link_libraries(model_tests PRIVATE igor::Model igor::Core Catch2::Catch2WithMain)
    target_compile_definitions(model_tests PRIVATE IGOR_MODELS_DIR="${CMAKE_SOURCE_DIR}/models")
    include(Catch)
    catch_discover_tests(model_tests)
endif()

# --- Core Module Tests (LegacyBridge, etc.) ---
file(GLOB CORE_MODULE_TEST_SOURCES "igor/Core/*.cpp")

if(CORE_MODULE_TEST_SOURCES)
    add_executable(core_module_tests ${CORE_MODULE_TEST_SOURCES})
    set_target_properties(core_module_tests PROPERTIES CXX_STANDARD 23)
    target_link_libraries(core_module_tests PRIVATE igor::Core igor::Model Catch2::Catch2WithMain)
    target_compile_definitions(core_module_tests PRIVATE IGOR_MODELS_DIR="${CMAKE_SOURCE_DIR}/models")
    include(Catch)
    catch_discover_tests(core_module_tests)
endif()
