enable_testing()
# CMake configuration for tests

# Find Catch2
find_package(Catch2 3 REQUIRED)

# On Windows, copy runtime DLLs (project shared libs + imported DLLs) next to
# test executables so the loader finds them in the application directory.
# This avoids 0xc0000142 (STATUS_DLL_NOT_FOUND) errors during Catch2 test
# discovery, which can occur when DL_PATHS / PATH manipulation is unreliable.
function(igor_copy_runtime_dlls target)
    if(WIN32)
        add_custom_command(TARGET ${target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_RUNTIME_DLLS:${target}> $<TARGET_FILE_DIR:${target}>
            COMMAND_EXPAND_LISTS
        )
    endif()
endfunction()

# Add test sources
file(GLOB TEST_SOURCES "*.cpp")

add_executable(igor_tests ${TEST_SOURCES})

# Link against the core library, model library and Catch2
target_link_libraries(igor_tests PRIVATE igor::Core igor::Model Catch2::Catch2WithMain)

# Provide the project source directory so tests can locate data files
target_compile_definitions(igor_tests PRIVATE IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

igor_copy_runtime_dlls(igor_tests)

# Set include directories if needed
# target_include_directories(igor_tests PRIVATE ../src/igor/Core)

# Register the tests with CTest
include(Catch)
catch_discover_tests(igor_tests
    PROPERTIES TIMEOUT 3600  # generation test with 1M sequences needs headroom
    DL_PATHS $<TARGET_FILE_DIR:igor_tests> $ENV{CONDA_PREFIX}/Library/bin
    DISCOVERY_MODE PRE_TEST
)

# --- Math Tests ---
file(GLOB MATH_TEST_SOURCES "igor/Math/*.cpp")

if(MATH_TEST_SOURCES)
    add_executable(math_tests ${MATH_TEST_SOURCES})
    set_target_properties(math_tests PROPERTIES CXX_STANDARD 23)
    target_link_libraries(math_tests PRIVATE igor::Math Catch2::Catch2WithMain)
    igor_copy_runtime_dlls(math_tests)
    include(Catch)
    catch_discover_tests(math_tests
        DL_PATHS $<TARGET_FILE_DIR:math_tests> $ENV{CONDA_PREFIX}/Library/bin
        DISCOVERY_MODE PRE_TEST
    )
endif()

# --- Model Tests ---
file(GLOB MODEL_TEST_SOURCES "igor/Model/*.cpp")

if(MODEL_TEST_SOURCES)
    add_executable(model_tests ${MODEL_TEST_SOURCES})
    set_target_properties(model_tests PROPERTIES CXX_STANDARD 23)
    target_link_libraries(model_tests PRIVATE igor::Model igor::Core Catch2::Catch2WithMain)
    target_compile_definitions(model_tests PRIVATE IGOR_MODELS_DIR="${CMAKE_SOURCE_DIR}/models")
    igor_copy_runtime_dlls(model_tests)
    include(Catch)
    catch_discover_tests(model_tests
        DL_PATHS $<TARGET_FILE_DIR:model_tests> $ENV{CONDA_PREFIX}/Library/bin
        DISCOVERY_MODE PRE_TEST
    )
endif()

# --- Core Module Tests (LegacyBridge, etc.) ---
file(GLOB CORE_MODULE_TEST_SOURCES "igor/Core/*.cpp")

if(CORE_MODULE_TEST_SOURCES)
    add_executable(core_module_tests ${CORE_MODULE_TEST_SOURCES})
    set_target_properties(core_module_tests PROPERTIES CXX_STANDARD 23)
    target_link_libraries(core_module_tests PRIVATE igor::Core igor::Model Catch2::Catch2WithMain)
    target_compile_definitions(core_module_tests PRIVATE IGOR_MODELS_DIR="${CMAKE_SOURCE_DIR}/models")
    igor_copy_runtime_dlls(core_module_tests)
    include(Catch)
    catch_discover_tests(core_module_tests
        DL_PATHS $<TARGET_FILE_DIR:core_module_tests> $ENV{CONDA_PREFIX}/Library/bin
        DISCOVERY_MODE PRE_TEST
    )
endif()
