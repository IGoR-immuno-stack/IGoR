### CMakeLists.txt ---

add_executable(core_tests
    test_dynamic_sequence_map.cpp
    test_layered_array.cpp
)

# Legacy tests: the tests of the code still under <Layer>/Legacy.
target_sources(core_tests
    PRIVATE
        Legacy/test_sanity.cpp
)

set_target_properties(core_tests PROPERTIES CXX_STANDARD 23)

target_link_libraries(core_tests
    PRIVATE
        igor::Core
        Catch2::Catch2WithMain)

target_compile_definitions(core_tests PRIVATE
    IGOR_SOURCE_DIR="${CMAKE_SOURCE_DIR}")

include(Catch)
# CI's main job runs `ctest -L unit`; the other two groups are run by the nightly job.
catch_discover_tests(core_tests
    TEST_SPEC "~[integration]~[slow]"
    TEST_PREFIX "unit::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "unit"
)

catch_discover_tests(core_tests
    TEST_SPEC "[integration]"
    TEST_PREFIX "integration::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "integration"
)

catch_discover_tests(core_tests
    TEST_SPEC "[convergence]"
    TEST_PREFIX "convergence::"
    PROPERTIES
        TIMEOUT 3600
        LABELS "convergence"
)

######################################################################
### CMakeLists.txt ends here
