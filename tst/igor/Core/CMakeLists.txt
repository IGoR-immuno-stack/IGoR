### CMakeLists.txt ---

set(CORE_TEST_SOURCES
    test_Contexts.cpp
    test_EventUtils.cpp
    test_fast_sampling.cpp
    test_generation.cpp
    test_inference.cpp
    test_sanity.cpp
    test_utils.cpp
)

add_executable(igor_tests ${CORE_TEST_SOURCES})

target_link_libraries(igor_tests
    PRIVATE
        igor::Core
        igor::Model
        Catch2::Catch2WithMain)

target_compile_definitions(igor_tests PRIVATE
    IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")


igor_copy_runtime_dlls(igor_tests)
include(Catch)
set(_igor_tests_dl_paths "$<TARGET_FILE_DIR:igor_tests>")
if(WIN32 AND DEFINED ENV{CONDA_PREFIX})
    list(APPEND _igor_tests_dl_paths "$ENV{CONDA_PREFIX}/Library/bin")
endif()

catch_discover_tests(igor_tests
    DISCOVERY_MODE PRE_TEST
    TEST_SPEC "~[integration]"
    TEST_PREFIX "unit::"
    DL_PATHS ${_igor_tests_dl_paths}
    PROPERTIES TIMEOUT 3600 LABELS "unit"
)

catch_discover_tests(igor_tests
    DISCOVERY_MODE PRE_TEST
    TEST_SPEC "[integration]"
    TEST_PREFIX "integration::"
    DL_PATHS ${_igor_tests_dl_paths}
    PROPERTIES TIMEOUT 3600 LABELS "integration"
)

######################################################################
### CMakeLists.txt ends here
