### CMakeLists.txt ---

set(CORE_TEST_SOURCES
    test_Contexts.cpp
    test_EventUtils.cpp
    test_fast_sampling.cpp
    test_generation.cpp
    test_inference.cpp
    test_poc.cpp
    test_sanity.cpp
    test_utils.cpp
)

add_executable(igor_tests ${CORE_TEST_SOURCES})

add_test(
    NAME igor_test_unit
    COMMAND igor_tests [unit] -g 0
)

add_test(
    NAME igor_test_integration
    COMMAND igor_tests [regression] -g 1
)

add_test(
    NAME igor_test_convergence
    COMMAND igor_tests [convergence] -g 2
)

add_test(
    NAME igor_test_regression
    COMMAND igor_tests [regression] -g 3
)

add_test(
    NAME igor_test_benchmark
    COMMAND igor_tests [benchmark] -g 4
)


target_link_libraries(igor_tests
    PRIVATE
        igor::Core
        Catch2::Catch2WithMain)

target_compile_definitions(igor_tests PRIVATE
    IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

include(Catch)
catch_discover_tests(igor_tests
    TEST_SPEC "~[integration]"
    TEST_PREFIX "unit::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "unit"
)

catch_discover_tests(igor_tests
    TEST_SPEC "[integration]"
    TEST_PREFIX "integration::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "integration"
)

######################################################################
### CMakeLists.txt ends here
