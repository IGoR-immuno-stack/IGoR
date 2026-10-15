### CMakeLists.txt ---

add_executable(alignment_tests)

# Legacy tests: the tests of the code still under <Layer>/Legacy.
target_sources(alignment_tests
    PRIVATE
        Legacy/test_aligner.cpp
        Legacy/test_aligner_benchmark.cpp
        Legacy/test_alignment_cigar.cpp
        Legacy/AlignerTestUtils.cpp
)

set_target_properties(alignment_tests PROPERTIES CXX_STANDARD 23)

target_link_libraries(alignment_tests
    PRIVATE
        igor::Alignment
        Catch2::Catch2WithMain)

target_compile_definitions(alignment_tests PRIVATE
    IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

# Whitebox tests reach into Alignment-internal Smith-Waterman DP symbols (see
# AlignerInternal.h) that aren't part of the layer's public API. ALIGNMENT_TESTING_ENABLED
# switches those symbols' ALIGNMENT_TESTING_EXPORT annotation to a real export/import (instead
# of hidden visibility) only for this test build, on both sides of the shared library boundary:
# Alignment itself (so it actually exports them) and alignment_tests (so it declares them the
# same way when importing).
target_compile_definitions(Alignment PRIVATE ALIGNMENT_TESTING_ENABLED)
target_compile_definitions(alignment_tests PRIVATE ALIGNMENT_TESTING_ENABLED)

include(Catch)
catch_discover_tests(alignment_tests
    TEST_SPEC "~[integration]~[slow]"
    TEST_PREFIX "unit::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "unit"
)

catch_discover_tests(alignment_tests
    TEST_SPEC "[integration]"
    TEST_PREFIX "integration::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "integration"
)

catch_discover_tests(alignment_tests
    TEST_SPEC "[convergence]"
    TEST_PREFIX "convergence::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "convergence"
)

######################################################################
### CMakeLists.txt ends here
