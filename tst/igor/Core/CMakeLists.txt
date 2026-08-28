### CMakeLists.txt ---

set(CORE_TEST_SOURCES
    test_aligner.cpp
    test_aligner_benchmark.cpp
    test_alignment_cigar.cpp
    test_Contexts.cpp
    test_EventUtils.cpp
    test_fast_sampling.cpp
    test_generation.cpp
    test_genetic_code.cpp
    test_inference.cpp
    test_memory_map.cpp
    test_model_format_v2.cpp
    test_sanity.cpp
    test_utils.cpp
    AlignerTestUtils.cpp
)

add_executable(igor_tests ${CORE_TEST_SOURCES})

target_link_libraries(igor_tests
    PRIVATE
        igor::Core
        Catch2::Catch2WithMain)

target_compile_definitions(igor_tests PRIVATE
    IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# Whitebox tests reach into Core-internal Smith-Waterman DP symbols (see AlignerInternal.h)
# that aren't part of Core's normal public API. CORE_TESTING_ENABLED switches those symbols'
# CORE_TESTING_EXPORT annotation to a real export/import (instead of hidden visibility) only
# for this test build, on both sides of the Core shared library boundary: Core itself (so it
# actually exports them) and igor_tests (so it declares them the same way when importing).
target_compile_definitions(Core PRIVATE CORE_TESTING_ENABLED)
target_compile_definitions(igor_tests PRIVATE CORE_TESTING_ENABLED)

include(Catch)
catch_discover_tests(igor_tests
    TEST_SPEC "~[integration]~[slow]"
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

catch_discover_tests(igor_tests
    TEST_SPEC "[convergence]"
    TEST_PREFIX "convergence::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "convergence"
)

######################################################################
### CMakeLists.txt ends here
